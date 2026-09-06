#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

/*
 * Independent full-frontier arithmetic validator for OEIS A255885.
 *
 * This is intentionally not a persistence-capable production program.  It
 * rebuilds its own SPF table, uses shared exact uint32_t counters with atomic
 * increments, obtains prime-power roots by Hensel lifting roots modulo p, and
 * tests omitted factors with the original exponent c-1.  The production
 * program instead uses per-thread saturated byte counters, primitive-root
 * lifting to p^(2e), and an order/valuation cascade.
 */

#define MAX_B 2000000000U
#define INDEX_CAP 254U
#define MAX_FACTORS 9U
#define MAX_TARGETS 230U
#define MAX_WITNESSES INDEX_CAP
#define TARGET_HASH_SIZE 512U
#define ROOT_CACHE_SIZE 16384U

__extension__ typedef unsigned __int128 u128;
__extension__ typedef __int128 i128;

typedef struct {
    uint32_t p;
    uint32_t e;
    uint32_t q;
    uint32_t h;
    uint64_t mod;
    uint64_t generator;
} LocalFactor;

typedef struct {
    uint64_t roots;
    uint64_t increments;
    uint64_t order_rejects;
    uint64_t direct_rejects;
    uint64_t composites;
} Metrics;

typedef struct {
    uint32_t cache_key[ROOT_CACHE_SIZE];
    uint32_t cache_value[ROOT_CACHE_SIZE];
    uint16_t witness_count[MAX_TARGETS];
    uint32_t *witness;
    Metrics metrics;
} ThreadState;

typedef struct {
    uint32_t n;
    uint32_t base;
} Target;

typedef struct {
    uint32_t key[TARGET_HASH_SIZE];
    uint16_t value[TARGET_HASH_SIZE];
} TargetMap;

typedef struct {
    uint32_t B;
    uint32_t c;
    uint32_t selected;
    unsigned nf;
    const LocalFactor *factor;
    uint32_t *counts;
    const TargetMap *targets;
    ThreadState *thread;
} EnumContext;

_Static_assert((TARGET_HASH_SIZE & (TARGET_HASH_SIZE - 1U)) == 0U,
               "target hash size must be a power of two");
_Static_assert((ROOT_CACHE_SIZE & (ROOT_CACHE_SIZE - 1U)) == 0U,
               "root cache size must be a power of two");
_Static_assert((uint64_t)MAX_B * (uint64_t)MAX_B <= UINT64_MAX,
               "B squared must fit uint64_t");

static _Noreturn void die(const char *message)
{
    (void)fprintf(stderr, "independent-validator fatal: %s\n", message);
    exit(EXIT_FAILURE);
}

static _Noreturn void die_errno(const char *message)
{
    (void)fprintf(stderr, "independent-validator fatal: %s: %s\n",
                  message, strerror(errno));
    exit(EXIT_FAILURE);
}

static double monotonic_seconds(void)
{
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) die_errno("clock_gettime");
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1.0e9;
}

static uint32_t gcd32(uint32_t a, uint32_t b)
{
    while (b != 0U) {
        uint32_t const r = a % b;
        a = b;
        b = r;
    }
    return a;
}

static uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)(((u128)a * (u128)b) % (u128)m);
}

static uint64_t pow_mod(uint64_t base, uint64_t exponent, uint64_t modulus)
{
    if (modulus < 2U) die("invalid modular exponentiation modulus");
    uint64_t result = 1U;
    base %= modulus;
    while (exponent != 0U) {
        if ((exponent & 1U) != 0U) result = mul_mod(result, base, modulus);
        exponent >>= 1U;
        if (exponent != 0U) base = mul_mod(base, base, modulus);
    }
    return result;
}

static uint64_t inverse_mod(uint64_t a, uint64_t modulus)
{
    if (modulus < 2U || a == 0U || a >= modulus)
        die("invalid modular inverse input");
    i128 old_t = 0;
    i128 t = 1;
    uint64_t old_r = modulus;
    uint64_t r = a;
    while (r != 0U) {
        uint64_t const q = old_r / r;
        uint64_t const nr = old_r - q * r;
        i128 const nt = old_t - (i128)q * t;
        old_r = r;
        r = nr;
        old_t = t;
        t = nt;
    }
    if (old_r != 1U) die("non-coprime modular inverse input");
    old_t %= (i128)modulus;
    if (old_t < 0) old_t += (i128)modulus;
    return (uint64_t)old_t;
}

static uint32_t *build_spf(uint32_t B)
{
    size_t const count = (size_t)B + 1U;
    if (count > SIZE_MAX / sizeof(uint32_t)) die("SPF allocation overflow");
    uint32_t *spf = calloc(count, sizeof(*spf));
    if (spf == NULL) die_errno("cannot allocate independent SPF table");
    if (B >= 1U) spf[1] = 1U;
    for (uint32_t p = 2U; p <= B / p; ++p) {
        if (spf[p] != 0U) continue;
        spf[p] = p;
        uint64_t const start = (uint64_t)p * (uint64_t)p;
        for (uint64_t j = start; j <= (uint64_t)B; j += (uint64_t)p) {
            if (spf[(size_t)j] == 0U) spf[(size_t)j] = p;
        }
    }
    for (uint32_t n = 2U; n <= B; ++n) {
        if (spf[n] == 0U) spf[n] = n;
    }
    return spf;
}

static unsigned factor_number(uint32_t n, const uint32_t *spf,
                              LocalFactor out[MAX_FACTORS])
{
    unsigned nf = 0U;
    while (n > 1U) {
        if (nf >= MAX_FACTORS) die("factor-count bound exceeded");
        uint32_t const p = spf[n];
        if (p < 2U || n % p != 0U) die("invalid independent SPF entry");
        uint32_t q = 1U;
        uint32_t e = 0U;
        do {
            if (q > UINT32_MAX / p) die("prime-power overflow");
            q *= p;
            n /= p;
            ++e;
        } while (n > 1U && spf[n] == p);
        out[nf].p = p;
        out[nf].q = q;
        out[nf].e = e;
        out[nf].h = 0U;
        out[nf].mod = (uint64_t)q * (uint64_t)q;
        out[nf].generator = 1U;
        ++nf;
    }
    return nf;
}

static unsigned distinct_factors(uint32_t n, const uint32_t *spf,
                                 uint32_t out[MAX_FACTORS])
{
    unsigned count = 0U;
    while (n > 1U) {
        if (count >= MAX_FACTORS) die("p-1 factor-count bound exceeded");
        uint32_t const p = spf[n];
        out[count++] = p;
        do n /= p; while (n > 1U && spf[n] == p);
    }
    return count;
}

static uint32_t primitive_root_mod_p(uint32_t p, const uint32_t *spf,
                                     ThreadState *state)
{
    uint32_t const slot = (uint32_t)(((uint64_t)p * UINT64_C(11400714819323198485))
                                      & (ROOT_CACHE_SIZE - 1U));
    if (state->cache_key[slot] == p) return state->cache_value[slot];
    uint32_t prime_divisors[MAX_FACTORS];
    unsigned const nd = distinct_factors(p - 1U, spf, prime_divisors);
    uint32_t g = 2U;
    for (;; ++g) {
        if (g >= p) die("primitive-root search failed");
        bool good = true;
        for (unsigned i = 0U; i < nd; ++i) {
            if (pow_mod(g, (uint64_t)(p - 1U) / prime_divisors[i], p) == 1U) {
                good = false;
                break;
            }
        }
        if (good) break;
    }
    state->cache_key[slot] = p;
    state->cache_value[slot] = g;
    return g;
}

static uint64_t hensel_subgroup_generator(uint32_t c, LocalFactor *factor,
                                          const uint32_t *spf,
                                          ThreadState *state)
{
    uint32_t const p = factor->p;
    uint32_t const h = factor->h;
    if (p == 2U || h == 1U) return 1U;

    uint32_t const primitive = primitive_root_mod_p(p, spf, state);
    uint64_t x = pow_mod(primitive, (p - 1U) / h, p);
    uint64_t modulus = p;
    unsigned const target_exponent = 2U * factor->e;
    uint64_t const k = (uint64_t)c - 1U;

    for (unsigned exponent = 1U; exponent < target_exponent; ++exponent) {
        if (modulus > UINT64_MAX / p) die("Hensel modulus overflow");
        uint64_t const next_modulus = modulus * p;
        uint64_t const fx_residue = pow_mod(x, k, next_modulus);
        if (fx_residue == 0U || (fx_residue - 1U) % modulus != 0U)
            die("invalid Hensel root invariant");
        uint64_t const quotient = (fx_residue - 1U) / modulus;
        uint64_t derivative = mul_mod(k % p, pow_mod(x % p, k - 1U, p), p);
        if (derivative == 0U) die("singular Hensel derivative");
        derivative = inverse_mod(derivative, p);
        uint64_t const correction = (p - mul_mod(quotient % p, derivative, p)) % p;
        x += correction * modulus;
        modulus = next_modulus;
    }
    if (modulus != factor->mod || pow_mod(x, k, modulus) != 1U)
        die("Hensel lift verification failed");
    if (pow_mod(x, h, modulus) != 1U) die("lifted subgroup order mismatch");
    return x;
}

static uint32_t choose_projection(uint32_t c, uint32_t B,
                                  const LocalFactor *factor, unsigned nf)
{
    uint32_t const full = (UINT32_C(1) << nf) - 1U;
    if ((uint64_t)c * (uint64_t)c <= B) return full;
    uint64_t best_cost = UINT64_MAX;
    uint64_t best_modulus = 0U;
    uint32_t best_mask = 0U;
    for (uint32_t mask = 1U; mask <= full; ++mask) {
        uint64_t d = 1U;
        uint64_t cost = 1U;
        for (unsigned i = 0U; i < nf; ++i) {
            if ((mask & (UINT32_C(1) << i)) == 0U) continue;
            d *= factor[i].q;
            cost *= factor[i].h;
        }
        if (d * d <= B) continue;
        /* Reverse tie-break versus the production implementation. */
        if (cost < best_cost ||
            (cost == best_cost && (d * d > best_modulus ||
                                   (d * d == best_modulus && mask > best_mask)))) {
            best_cost = cost;
            best_modulus = d * d;
            best_mask = mask;
        }
    }
    if (best_mask == 0U) die("projection subset search failed");
    return best_mask;
}

static int target_lookup(const TargetMap *map, uint32_t base)
{
    uint32_t slot = (uint32_t)(((uint64_t)base * UINT64_C(2654435761))
                                & (TARGET_HASH_SIZE - 1U));
    for (uint32_t probes = 0U; probes < TARGET_HASH_SIZE; ++probes) {
        if (map->key[slot] == 0U) return -1;
        if (map->key[slot] == base) return (int)map->value[slot] - 1;
        slot = (slot + 1U) & (TARGET_HASH_SIZE - 1U);
    }
    die("target hash lookup exhausted");
}

static void record_incidence(EnumContext *ctx, uint32_t base)
{
#pragma omp atomic update
    ctx->counts[base] += 1U;
    ++ctx->thread->metrics.increments;

    int const target = target_lookup(ctx->targets, base);
    if (target >= 0) {
        uint16_t *used = &ctx->thread->witness_count[(unsigned)target];
        if (*used >= MAX_WITNESSES) die("thread-local witness capacity exceeded");
        size_t const offset = (size_t)(unsigned)target * MAX_WITNESSES + *used;
        ctx->thread->witness[offset] = ctx->c;
        ++*used;
    }
}

static bool omitted_factors_hold(uint64_t base, uint32_t c,
                                 const LocalFactor *factor, unsigned nf,
                                 uint32_t selected, Metrics *metrics)
{
    uint64_t const k = (uint64_t)c - 1U;
    for (unsigned i = 0U; i < nf; ++i) {
        if ((selected & (UINT32_C(1) << i)) != 0U) continue;
        uint32_t const p = factor[i].p;
        if (p != 2U) {
            uint32_t const h = gcd32(c - 1U, p - 1U);
            if (pow_mod(base, h, p) != 1U) {
                ++metrics->order_rejects;
                return false;
            }
        }
        if (pow_mod(base, k, factor[i].mod) != 1U) {
            ++metrics->direct_rejects;
            return false;
        }
    }
    return true;
}

static void accept_root(EnumContext *ctx, uint64_t root, uint64_t modulus)
{
    ++ctx->thread->metrics.roots;
    if (modulus <= ctx->B) {
        uint64_t base = root;
        if (base <= ctx->c) {
            uint64_t const steps = ((uint64_t)ctx->c - base) / modulus + 1U;
            base += steps * modulus;
        }
        while (base <= ctx->B) {
            record_incidence(ctx, (uint32_t)base);
            if (modulus > (uint64_t)ctx->B - base) break;
            base += modulus;
        }
    } else if (root > ctx->c && root <= ctx->B &&
               omitted_factors_hold(root, ctx->c, ctx->factor, ctx->nf,
                                    ctx->selected, &ctx->thread->metrics)) {
        record_incidence(ctx, (uint32_t)root);
    }
}

static void enumerate_roots(EnumContext *ctx, unsigned start,
                            uint64_t residue, uint64_t modulus)
{
    unsigned i = start;
    while (i < ctx->nf && (ctx->selected & (UINT32_C(1) << i)) == 0U) ++i;
    if (i == ctx->nf) {
        accept_root(ctx, residue, modulus);
        return;
    }

    LocalFactor const *local = &ctx->factor[i];
    uint64_t root = 1U;
    uint64_t const reduced = modulus == 1U ? 0U : residue % local->mod;
    uint64_t const inverse = modulus == 1U ? 0U :
        inverse_mod(modulus % local->mod, local->mod);
    for (uint32_t j = 0U; j < local->h; ++j) {
        uint64_t combined;
        uint64_t combined_modulus;
        if (modulus == 1U) {
            combined = root;
            combined_modulus = local->mod;
        } else {
            uint64_t const difference = root >= reduced ? root - reduced :
                local->mod - (reduced - root);
            uint64_t const t = mul_mod(difference, inverse, local->mod);
            u128 const wide = (u128)residue + (u128)modulus * t;
            u128 const wide_modulus = (u128)modulus * local->mod;
            if (wide_modulus > UINT64_MAX) die("CRT modulus overflow");
            combined_modulus = (uint64_t)wide_modulus;
            combined = (uint64_t)(wide % wide_modulus);
        }
        enumerate_roots(ctx, i + 1U, combined, combined_modulus);
        root = mul_mod(root, local->generator, local->mod);
    }
    if (root != 1U) die("local root cycle did not close");
}

static void process_composite(uint32_t c, uint32_t B, const uint32_t *spf,
                              uint32_t *counts, const TargetMap *targets,
                              ThreadState *thread)
{
    LocalFactor factor[MAX_FACTORS];
    unsigned const nf = factor_number(c, spf, factor);
    if (nf == 1U && factor[0].e == 1U) die("prime reached composite path");
    for (unsigned i = 0U; i < nf; ++i) {
        factor[i].h = factor[i].p == 2U ? 1U : gcd32(c - 1U, factor[i].p - 1U);
    }
    uint32_t const selected = choose_projection(c, B, factor, nf);
    for (unsigned i = 0U; i < nf; ++i) {
        if ((selected & (UINT32_C(1) << i)) != 0U)
            factor[i].generator = hensel_subgroup_generator(c, &factor[i], spf, thread);
    }
    ++thread->metrics.composites;
    EnumContext ctx = {B, c, selected, nf, factor, counts, targets, thread};
    enumerate_roots(&ctx, 0U, 0U, 1U);
}

static uint32_t parse_u32(const char *text, const char *field)
{
    if (text[0] < '0' || text[0] > '9') {
        (void)fprintf(stderr, "invalid %s\n", field);
        exit(EXIT_FAILURE);
    }
    for (const char *p = text; *p != '\0' && *p != '\n'; ++p) {
        if (*p < '0' || *p > '9') {
            (void)fprintf(stderr, "invalid %s\n", field);
            exit(EXIT_FAILURE);
        }
    }
    errno = 0;
    char *end = NULL;
    unsigned long long const value = strtoull(text, &end, 10);
    if (errno != 0 || end == text ||
        (*end != '\0' && !(*end == '\n' && end[1] == '\0')) ||
        value > UINT32_MAX) {
        (void)fprintf(stderr, "invalid %s\n", field);
        exit(EXIT_FAILURE);
    }
    return (uint32_t)value;
}

static uint32_t read_results(const char *path, uint32_t expected[INDEX_CAP + 1U],
                             Target targets[MAX_TARGETS], unsigned *target_count)
{
    FILE *file = fopen(path, "r");
    if (file == NULL) die_errno("cannot open results.tsv");
    char *line = NULL;
    size_t capacity = 0U;
    ssize_t length;
    uint32_t B = 0U;
    unsigned rows = 0U;
    *target_count = 0U;
    for (unsigned i = 0U; i <= INDEX_CAP; ++i) expected[i] = UINT32_MAX;
    while ((length = getline(&line, &capacity, file)) >= 0) {
        (void)length;
        if (strncmp(line, "B\t", 2U) == 0) {
            B = parse_u32(line + 2U, "results B");
            continue;
        }
        if (line[0] < '0' || line[0] > '9') continue;
        char *tab = strchr(line, '\t');
        if (tab == NULL) die("malformed numeric results row");
        *tab = '\0';
        uint32_t const n = parse_u32(line, "result index");
        char *base_text = tab + 1U;
        char *newline = strchr(base_text, '\n');
        if (newline != NULL) *newline = '\0';
        if (n != rows + 1U || n > INDEX_CAP) die("result row order mismatch");
        if (strcmp(base_text, "NA") != 0) {
            uint32_t const base = parse_u32(base_text, "result base");
            if (base < 2U || base > B) die("result base outside bound");
            expected[n] = base;
            if (n >= 25U) {
                if (*target_count >= MAX_TARGETS) die("too many validation targets");
                targets[*target_count].n = n;
                targets[*target_count].base = base;
                ++*target_count;
            }
        }
        ++rows;
    }
    free(line);
    if (fclose(file) != 0) die_errno("cannot close results.tsv");
    if (B < 4U || B > MAX_B || rows != INDEX_CAP)
        die("results metadata/domain mismatch");
    return B;
}

static void build_target_map(const Target *targets, unsigned count, TargetMap *map)
{
    memset(map, 0, sizeof(*map));
    for (unsigned i = 0U; i < count; ++i) {
        uint32_t const base = targets[i].base;
        uint32_t slot = (uint32_t)(((uint64_t)base * UINT64_C(2654435761))
                                   & (TARGET_HASH_SIZE - 1U));
        while (map->key[slot] != 0U) {
            if (map->key[slot] == base) die("duplicate candidate base");
            slot = (slot + 1U) & (TARGET_HASH_SIZE - 1U);
        }
        map->key[slot] = base;
        map->value[slot] = (uint16_t)(i + 1U);
    }
}

static Metrics recompute(uint32_t B, const uint32_t *spf, uint32_t *counts,
                         const TargetMap *targets, ThreadState *states,
                         unsigned requested_threads, unsigned *actual_threads)
{
    double const start = monotonic_seconds();
    *actual_threads = 0U;
#pragma omp parallel num_threads(requested_threads) default(none) \
    shared(B, spf, counts, targets, states, actual_threads, start, stderr)
    {
        unsigned const tid = (unsigned)omp_get_thread_num();
#pragma omp single
        *actual_threads = (unsigned)omp_get_num_threads();
#pragma omp for schedule(dynamic, 1024)
        for (uint64_t wide_c = 4U; wide_c < (uint64_t)B; ++wide_c) {
            uint32_t const c = (uint32_t)wide_c;
            if (spf[c] != c)
                process_composite(c, B, spf, counts, targets, &states[tid]);
            if (c % 50000000U == 0U) {
#pragma omp critical(a255885_independent_progress)
                {
                    double const elapsed = monotonic_seconds() - start;
                    (void)fprintf(stderr,
                        "independent-heartbeat\tc=%u\tfraction=%.6f\telapsed=%.3f\n",
                        c, (double)c / (double)B, elapsed);
                }
            }
        }
    }
    Metrics total = {0U, 0U, 0U, 0U, 0U};
    for (unsigned t = 0U; t < *actual_threads; ++t) {
        total.roots += states[t].metrics.roots;
        total.increments += states[t].metrics.increments;
        total.order_rejects += states[t].metrics.order_rejects;
        total.direct_rejects += states[t].metrics.direct_rejects;
        total.composites += states[t].metrics.composites;
    }
    return total;
}

static uint32_t load_le32(const unsigned char *p)
{
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8U) |
           ((uint32_t)p[2] << 16U) | ((uint32_t)p[3] << 24U);
}

static uint64_t load_le64(const unsigned char *p)
{
    uint64_t value = 0U;
    for (unsigned i = 0U; i < 8U; ++i) value |= (uint64_t)p[i] << (8U * i);
    return value;
}

static void compare_checkpoint(const char *path, uint32_t B, const uint32_t *counts)
{
    int const fd = open(path, O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
    if (fd < 0) die_errno("cannot open frozen count checkpoint");
    struct stat st;
    if (fstat(fd, &st) != 0) die_errno("cannot stat frozen count checkpoint");
    uint64_t const expected_size = 64U + (uint64_t)B + 1U;
    if (st.st_size < 0 || (uint64_t)st.st_size != expected_size)
        die("frozen count checkpoint size mismatch");
    void *mapping = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (mapping == MAP_FAILED) die_errno("cannot map frozen count checkpoint");
    unsigned char const *data = mapping;
    if (memcmp(data, "A255885-BINARY", 15U) != 0 ||
        load_le32(data + 16U) != 1U || load_le32(data + 20U) != 2U ||
        load_le32(data + 24U) != B || load_le32(data + 28U) != 1U ||
        load_le64(data + 32U) != (uint64_t)B + 1U)
        die("frozen count checkpoint header mismatch");
    unsigned char const *payload = data + 64U;
    for (uint32_t b = 0U; b <= B; ++b) {
        if (counts[b] > UINT8_MAX || payload[b] != (unsigned char)counts[b]) {
            (void)fprintf(stderr,
                "counter mismatch at b=%u: independent=%u frozen=%u\n",
                b, counts[b], (unsigned)payload[b]);
            exit(EXIT_FAILURE);
        }
    }
    if (munmap(mapping, (size_t)st.st_size) != 0) die_errno("cannot unmap checkpoint");
    if (close(fd) != 0) die_errno("cannot close checkpoint");
}

static int compare_u32(const void *a, const void *b)
{
    uint32_t const av = *(const uint32_t *)a;
    uint32_t const bv = *(const uint32_t *)b;
    return av < bv ? -1 : av > bv ? 1 : 0;
}

static void write_certificate(const char *path, uint32_t B,
                              const uint32_t expected[INDEX_CAP + 1U],
                              const Target *targets, unsigned target_count,
                              const uint32_t minima[INDEX_CAP + 1U],
                              const ThreadState *states, unsigned threads,
                              Metrics metrics, double elapsed)
{
    if (access(path, F_OK) == 0) die("certificate target already exists");
    size_t const length = strlen(path);
    if (length > SIZE_MAX - 5U) die("certificate path too long");
    char *temporary = malloc(length + 5U);
    if (temporary == NULL) die_errno("cannot allocate certificate path");
    (void)memcpy(temporary, path, length);
    (void)memcpy(temporary + length, ".tmp", 5U);
    int const fd = open(temporary, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0644);
    if (fd < 0) die_errno("cannot create certificate temporary file");
    FILE *file = fdopen(fd, "w");
    if (file == NULL) die_errno("cannot open certificate stream");

    (void)fprintf(file,
        "A255885-INDEPENDENT-CERT\t1\n"
        "B\t%u\nthreads\t%u\nroots\t%" PRIu64 "\nincrements\t%" PRIu64
        "\norder_rejects\t%" PRIu64 "\ndirect_rejects\t%" PRIu64
        "\ncomposites\t%" PRIu64 "\nelapsed_seconds\t%.9f\n"
        "full_counter_comparison\tPASS\nminima_comparison\tPASS\n"
        "new_target_count\t%u\n",
        B, threads, metrics.roots, metrics.increments, metrics.order_rejects,
        metrics.direct_rejects, metrics.composites, elapsed, target_count);

    for (uint32_t n = 1U; n <= INDEX_CAP; ++n) {
        if (expected[n] == UINT32_MAX) {
            (void)fprintf(file, "M\t%u\tNA\n", n);
        } else {
            (void)fprintf(file, "M\t%u\t%u\n", n, minima[n]);
        }
    }

    uint32_t merged[MAX_WITNESSES];
    for (unsigned slot = 0U; slot < target_count; ++slot) {
        unsigned used = 0U;
        for (unsigned t = 0U; t < threads; ++t) {
            unsigned const local_count = states[t].witness_count[slot];
            if (used + local_count > MAX_WITNESSES) die("merged witness capacity exceeded");
            size_t const offset = (size_t)slot * MAX_WITNESSES;
            for (unsigned j = 0U; j < local_count; ++j)
                merged[used++] = states[t].witness[offset + j];
        }
        qsort(merged, used, sizeof(*merged), compare_u32);
        if (used != targets[slot].n) {
            (void)fprintf(stderr,
                "witness cardinality mismatch for n=%u: got %u\n",
                targets[slot].n, used);
            exit(EXIT_FAILURE);
        }
        for (unsigned j = 0U; j < used; ++j) {
            if ((j > 0U && merged[j] == merged[j - 1U]) ||
                merged[j] >= targets[slot].base ||
                pow_mod(targets[slot].base, (uint64_t)merged[j] - 1U,
                        (uint64_t)merged[j] * merged[j]) != 1U)
                die("direct witness verification failed");
        }
        (void)fprintf(file, "W\t%u\t%u\t%u\t",
                      targets[slot].n, targets[slot].base, used);
        for (unsigned j = 0U; j < used; ++j)
            (void)fprintf(file, "%s%u", j == 0U ? "" : ",", merged[j]);
        (void)fputc('\n', file);
    }

    if (fflush(file) != 0 || fsync(fd) != 0) die_errno("cannot flush certificate");
    if (fclose(file) != 0) die_errno("cannot close certificate");
    if (rename(temporary, path) != 0) die_errno("cannot commit certificate");
    free(temporary);
}

static void initialize_states(ThreadState *states, unsigned threads)
{
    for (unsigned t = 0U; t < threads; ++t) {
        size_t const elements = (size_t)MAX_TARGETS * MAX_WITNESSES;
        states[t].witness = calloc(elements, sizeof(uint32_t));
        if (states[t].witness == NULL) die_errno("cannot allocate witness storage");
    }
}

static void free_states(ThreadState *states, unsigned threads)
{
    for (unsigned t = 0U; t < threads; ++t) free(states[t].witness);
}

static void run_self_test(uint32_t B)
{
    if (B < 20U || B > 5000U) die("self-test bound must be in [20,5000]");
    uint32_t *spf = build_spf(B);
    uint32_t *counts = calloc((size_t)B + 1U, sizeof(*counts));
    uint32_t *direct = calloc((size_t)B + 1U, sizeof(*direct));
    if (counts == NULL || direct == NULL) die_errno("self-test allocation failed");
    TargetMap empty;
    memset(&empty, 0, sizeof(empty));
    unsigned const requested = 4U;
    ThreadState *states = calloc(requested, sizeof(*states));
    if (states == NULL) die_errno("self-test state allocation failed");
    initialize_states(states, requested);
    unsigned actual = 0U;
    (void)recompute(B, spf, counts, &empty, states, requested, &actual);

    for (uint32_t c = 4U; c < B; ++c) {
        if (spf[c] == c) continue;
        uint64_t const modulus = (uint64_t)c * c;
        for (uint32_t b = c + 1U; b <= B; ++b) {
            if (pow_mod(b, (uint64_t)c - 1U, modulus) == 1U) ++direct[b];
        }
    }
    for (uint32_t b = 0U; b <= B; ++b) {
        if (counts[b] != direct[b]) {
            (void)fprintf(stderr, "self-test mismatch at b=%u: %u != %u\n",
                          b, counts[b], direct[b]);
            exit(EXIT_FAILURE);
        }
    }
    free_states(states, requested);
    free(states);
    free(direct);
    free(counts);
    free(spf);
    (void)printf("INDEPENDENT SELF-TEST PASS\tB=%u\tthreads=%u\n", B, actual);
}

static void run_validation(const char *results_path, const char *counts_path,
                           const char *certificate_path, unsigned requested_threads)
{
    uint32_t expected[INDEX_CAP + 1U];
    Target targets[MAX_TARGETS];
    unsigned target_count = 0U;
    uint32_t const B = read_results(results_path, expected, targets, &target_count);
    TargetMap target_map;
    build_target_map(targets, target_count, &target_map);

    (void)fprintf(stderr, "independent-start\tB=%u\ttargets=%u\trequested_threads=%u\n",
                  B, target_count, requested_threads);
    double const total_start = monotonic_seconds();
    uint32_t *spf = build_spf(B);
    (void)fprintf(stderr, "independent-stage\tspf=complete\telapsed=%.3f\n",
                  monotonic_seconds() - total_start);
    size_t const count = (size_t)B + 1U;
    uint32_t *counts = calloc(count, sizeof(*counts));
    if (counts == NULL) die_errno("cannot allocate exact counter array");
    ThreadState *states = calloc(requested_threads, sizeof(*states));
    if (states == NULL) die_errno("cannot allocate thread states");
    initialize_states(states, requested_threads);

    unsigned actual_threads = 0U;
    double const arithmetic_start = monotonic_seconds();
    Metrics const metrics = recompute(B, spf, counts, &target_map, states,
                                      requested_threads, &actual_threads);
    double const arithmetic_elapsed = monotonic_seconds() - arithmetic_start;
    (void)fprintf(stderr,
        "independent-stage\tarithmetic=complete\tthreads=%u\telapsed=%.3f\n",
        actual_threads, arithmetic_elapsed);

    compare_checkpoint(counts_path, B, counts);
    (void)fprintf(stderr, "independent-stage\tfull_counter_comparison=PASS\n");

    uint32_t minima[INDEX_CAP + 1U];
    memset(minima, 0, sizeof(minima));
    for (uint32_t b = 2U; b <= B; ++b) {
        uint32_t const value = counts[b];
        if (value >= 1U && value <= INDEX_CAP && minima[value] == 0U)
            minima[value] = b;
    }
    for (uint32_t n = 1U; n <= INDEX_CAP; ++n) {
        uint32_t const independent = minima[n] == 0U ? UINT32_MAX : minima[n];
        if (independent != expected[n]) {
            (void)fprintf(stderr,
                "minimum mismatch at n=%u: independent=%s%u frozen=%s%u\n",
                n, minima[n] == 0U ? "NA" : "", minima[n],
                expected[n] == UINT32_MAX ? "NA" : "", expected[n]);
            exit(EXIT_FAILURE);
        }
    }
    write_certificate(certificate_path, B, expected, targets, target_count,
                      minima, states, actual_threads, metrics, arithmetic_elapsed);
    (void)fprintf(stderr, "independent-stage\tcertificate=committed\telapsed=%.3f\n",
                  monotonic_seconds() - total_start);
    free_states(states, requested_threads);
    free(states);
    free(counts);
    free(spf);
    (void)printf("INDEPENDENT FULL VALIDATION PASS\tB=%u\ttargets=%u\tthreads=%u\n",
                 B, target_count, actual_threads);
}

static unsigned parse_threads(const char *text)
{
    uint32_t const value = parse_u32(text, "thread count");
    if (value < 1U || value > 64U) die("thread count must be in [1,64]");
    return (unsigned)value;
}

int main(int argc, char **argv)
{
    if (argc == 3 && strcmp(argv[1], "self-test") == 0) {
        run_self_test(parse_u32(argv[2], "self-test bound"));
        return EXIT_SUCCESS;
    }
    if (argc == 7 && strcmp(argv[1], "validate") == 0 &&
        strcmp(argv[5], "--threads") == 0) {
        run_validation(argv[2], argv[3], argv[4], parse_threads(argv[6]));
        return EXIT_SUCCESS;
    }
    (void)fprintf(stderr,
        "Usage:\n"
        "  a255885-independent self-test B\n"
        "  a255885-independent validate RESULTS_TSV ACTIVE_COUNT_BIN CERTIFICATE --threads N\n");
    return EXIT_FAILURE;
}
