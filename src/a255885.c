/*
 * A255885 - interval-aware inverse root sieve (external v6 candidate)
 *
 * Copyright (c) 2026 Carlo Corti.
 * SPDX-License-Identifier: MIT
 *
 * For 2 <= b <= B this program computes
 *
 *   C(b) = #{ composite c : 2 <= c < b and b^(c-1) == 1 (mod c^2) }
 *
 * and, only after complete coverage, extracts the least b with C(b)=n for
 * 1 <= n <= 254.  Counters saturate at 255, which preserves all equalities
 * C(b)=n in the supported target domain.
 *
 * Prior art: Robert Israel's program for OEIS A273785 already uses the
 * general inversion "composite modulus -> roots -> base progressions".
 * This implementation adds an exact interval-aware complete-prime-power
 * divisor projection, bounded C17/OpenMP engineering, and fail-closed
 * segmented persistence for A255885.
 *
 * Mathematical kernel.  Write c=prod p_i^e_i and k=c-1.  CRT reduces the
 * congruence modulo c^2 to the coprime moduli p_i^(2e_i).  Since p_i does
 * not divide k, the local kernel has size 1 at p=2 and
 * gcd(k,p-1) at odd p.  For odd p, a primitive root modulo p^2 is also a
 * primitive root modulo every higher p-power; raising it to phi/h generates
 * the h-element local kernel.  For p=2, exponentiation by odd k is an
 * automorphism of the 2-group of units, so the only local root is 1.
 *
 * If c^2>B, a divisor d made of complete prime-power factors is selected
 * with d^2>B and minimum projected kernel size (numeric factor-mask tie
 * break).  Since b<=B<d^2, an admitted full solution equals its projected
 * residue.  Every omitted factor is then checked exactly.  For odd p^e,
 *
 *   b^k == 1 (mod p^(2e))
 *   iff b^gcd(k,p-1) == 1 (mod p)
 *       and b^(p-1) == 1 (mod p^(2e)).
 *
 * The order condition is the first test.  The second is the exact p-adic
 * valuation condition (Agoh-Dilcher-Skula, Lemma 5.1 / Corollary 5.2).
 * For p=2 the exact test is b == 1 (mod 2^(2e)).
 *
 * Persistence is Linux/POSIX specific by design.  Scientific arithmetic is
 * C17 plus OpenMP; no CAS or external executable is used.  The verify command
 * certifies identity, persistence, coverage, and result derivation; it does
 * not independently recompute the arithmetic definition.
 */

#define _POSIX_C_SOURCE 200809L

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <signal.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#ifndef SOURCE_SHA256
#define SOURCE_SHA256 "UNSET"
#endif
#ifndef MAKEFILE_SHA256
#define MAKEFILE_SHA256 "UNSET"
#endif

#define PROGRAM_ID "A255885"
#define PROGRAM_VERSION 1U
#define CONFIG_MAGIC "A255885-CONFIG"
#define MANIFEST_MAGIC "A255885-MANIFEST"
#define BINARY_MAGIC "A255885-BINARY"
#define FORMAT_VERSION 1U
#define MAX_B 2000000000U
#define INDEX_CAP 254U
#define COUNTER_SAT 255U
#define MAX_FACTORS 9U
#define MAX_PATH_LEN 4096U
#define HASH_HEX_LEN 64U
#define BIN_HEADER_SIZE 64U
#define RESULTS_BUFFER_SIZE 8192U
#define ROOT_CACHE_SIZE 8192U
#define ARTIFACT_COUNT 7U
#define CACHE_LINE_SIZE 64U
#define OMP_DYNAMIC_CHUNK 64U
#define HEARTBEAT_INTERVAL_NS 5000000000ULL
#define HEARTBEAT_ROOT_STRIDE 65536U
#define HEARTBEAT_COMPOSITE_STRIDE 1024U

static bool g_test_fault_armed=false;

__extension__ typedef unsigned __int128 u128;
__extension__ typedef __int128 i128;

_Static_assert(COUNTER_SAT==INDEX_CAP+1U,
               "counter saturation must be the first unsupported index");
_Static_assert(COUNTER_SAT<=UINT8_MAX,
               "counter saturation must fit the checkpoint element type");
_Static_assert(MAX_FACTORS==9U && 223092870ULL<=(uint64_t)MAX_B &&
               6469693230ULL>(uint64_t)MAX_B,
               "factor capacity must cover every admitted integer");
_Static_assert(MAX_FACTORS<sizeof(unsigned)*CHAR_BIT,
               "factor masks must fit unsigned");
_Static_assert(MAX_B<=UINT32_MAX-16U,
               "B and default-segment next values must not wrap uint32_t");
_Static_assert((uint64_t)MAX_B*(uint64_t)MAX_B<=UINT64_MAX,
               "B squared must fit uint64_t");
_Static_assert(256U+16U*INDEX_CAP<=RESULTS_BUFFER_SIZE,
               "results buffer must cover every supported index row");

typedef struct {
    uint32_t p;
    uint32_t e;
    uint32_t q;
    uint32_t h;
    uint64_t modulus;
    uint64_t generator;
} Factor;

typedef struct {
    uint64_t roots_generated;
    uint64_t range_rejects;
    uint64_t order_rejects;
    uint64_t valuation_tests;
    uint64_t valuation_rejects;
    uint64_t two_adic_tests;
    uint64_t two_adic_rejects;
    uint64_t increments;
    uint64_t saturations;
    uint64_t composites;
} Metrics;

typedef struct {
    _Alignas(CACHE_LINE_SIZE) Metrics value;
} MetricsSlot;

typedef struct {
    uint32_t key[ROOT_CACHE_SIZE];
    uint64_t value[ROOT_CACHE_SIZE];
} RootCache;

typedef struct {
    _Alignas(CACHE_LINE_SIZE) RootCache value;
} RootCacheSlot;

_Static_assert(sizeof(MetricsSlot)%CACHE_LINE_SIZE==0U,
               "Metrics slots must be cache-line sized");
_Static_assert(_Alignof(MetricsSlot)>=CACHE_LINE_SIZE,
               "Metrics slots must be cache-line aligned");
_Static_assert(sizeof(RootCacheSlot)%CACHE_LINE_SIZE==0U,
               "Root-cache slots must be cache-line sized");
_Static_assert(_Alignof(RootCacheSlot)>=CACHE_LINE_SIZE,
               "Root-cache slots must be cache-line aligned");
_Static_assert((ROOT_CACHE_SIZE&(ROOT_CACHE_SIZE-1U))==0U,
               "Root cache size must be a power of two");

typedef struct {
    uint32_t B;
    uint32_t segment_span;
    uint32_t test_noncampaign;
    char source_hash[HASH_HEX_LEN + 1U];
    char makefile_hash[HASH_HEX_LEN + 1U];
    char binary_hash[HASH_HEX_LEN + 1U];
} Config;

typedef struct {
    uint64_t segment;
    uint32_t c_lo;
    uint32_t c_hi;
    uint32_t old_frontier;
    uint32_t frontier;
    unsigned threads;
    unsigned active_slot;
    char state[16];
    char config_hash[HASH_HEX_LEN + 1U];
    char spf_hash[HASH_HEX_LEN + 1U];
    char count_hash[HASH_HEX_LEN + 1U];
    Metrics metrics;
    double elapsed;
    char stop_reason[32];
} ManifestRow;

typedef struct {
    atomic_uint_fast64_t completed;
    atomic_uint_fast64_t next_ns;
    uint64_t start_ns;
    uint64_t total;
    uint32_t c_lo;
    uint32_t c_hi;
} Heartbeat;

_Noreturn static void die(const char *fmt, ...) __attribute__((format(printf, 1, 2)));

_Noreturn static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    (void)fprintf(stderr, "fatal: ");
    (void)vfprintf(stderr, fmt, ap);
    (void)fprintf(stderr, "\n");
    va_end(ap);
    exit(EXIT_FAILURE);
}

static void test_fault_exit(const char *stage, int code)
{
    const char *requested=getenv("A255885_TEST_FAULT");
    if(g_test_fault_armed && requested!=NULL && strcmp(requested,stage)==0) {
        (void)fprintf(stderr,"test_fault\tstage=%s\n",stage);
        _exit(code);
    }
}

static bool add_u64(uint64_t a, uint64_t b, uint64_t *out)
{
    if (UINT64_MAX - a < b) {
        return false;
    }
    *out = a + b;
    return true;
}

static bool mul_size(size_t a, size_t b, size_t *out)
{
    if (a != 0U && b > SIZE_MAX / a) {
        return false;
    }
    *out = a * b;
    return true;
}

static uint32_t gcd_u32(uint32_t a, uint32_t b)
{
    while (b != 0U) {
        uint32_t const t = a % b;
        a = b;
        b = t;
    }
    return a;
}

static uint64_t mul_mod_u64(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)(((u128)a * (u128)b) % (u128)m);
}

static uint64_t pow_mod_u64(uint64_t base, uint64_t exp, uint64_t mod)
{
    if (mod==0U) die("zero modular-arithmetic modulus");
    uint64_t result = 1U % mod;
    base %= mod;
    while (exp != 0U) {
        if ((exp & 1U) != 0U) {
            result = mul_mod_u64(result, base, mod);
        }
        exp >>= 1U;
        if (exp != 0U) {
            base = mul_mod_u64(base, base, mod);
        }
    }
    return result;
}

static uint64_t inv_mod_u64(uint64_t a, uint64_t m)
{
    if (m<2U) die("invalid inverse modulus");
    i128 old_r = (i128)a;
    i128 r = (i128)m;
    i128 old_s = 1;
    i128 s = 0;
    while (r != 0) {
        i128 const q = old_r / r;
        i128 const nr = old_r - q * r;
        i128 const ns = old_s - q * s;
        old_r = r;
        r = nr;
        old_s = s;
        s = ns;
    }
    if (old_r != 1) {
        die("internal CRT inverse does not exist");
    }
    old_s %= (i128)m;
    if (old_s < 0) {
        old_s += (i128)m;
    }
    return (uint64_t)old_s;
}

/* Minimal SHA-256 implementation used for all persisted identities. */
typedef struct {
    uint32_t state[8];
    uint64_t bitlen;
    uint8_t data[64];
    size_t datalen;
} Sha256;

static uint32_t rotr32(uint32_t x, unsigned n)
{
    return (x >> n) | (x << (32U - n));
}

static void sha256_transform(Sha256 *ctx, const uint8_t block[64])
{
    static const uint32_t k[64] = {
        0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,
        0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
        0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,
        0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
        0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,
        0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
        0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,
        0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
        0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,
        0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
        0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,
        0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
        0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,
        0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
        0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,
        0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U
    };
    uint32_t w[64];
    for (unsigned i = 0U; i < 16U; ++i) {
        size_t const j = (size_t)i * 4U;
        w[i] = ((uint32_t)block[j] << 24U) |
               ((uint32_t)block[j + 1U] << 16U) |
               ((uint32_t)block[j + 2U] << 8U) |
               (uint32_t)block[j + 3U];
    }
    for (unsigned i = 16U; i < 64U; ++i) {
        uint32_t const s0 = rotr32(w[i-15U],7U) ^ rotr32(w[i-15U],18U) ^ (w[i-15U] >> 3U);
        uint32_t const s1 = rotr32(w[i-2U],17U) ^ rotr32(w[i-2U],19U) ^ (w[i-2U] >> 10U);
        w[i] = w[i-16U] + s0 + w[i-7U] + s1;
    }
    uint32_t a=ctx->state[0], b=ctx->state[1], c=ctx->state[2], d=ctx->state[3];
    uint32_t e=ctx->state[4], f=ctx->state[5], g=ctx->state[6], h=ctx->state[7];
    for (unsigned i = 0U; i < 64U; ++i) {
        uint32_t const s1 = rotr32(e,6U) ^ rotr32(e,11U) ^ rotr32(e,25U);
        uint32_t const ch = (e & f) ^ ((~e) & g);
        uint32_t const t1 = h + s1 + ch + k[i] + w[i];
        uint32_t const s0 = rotr32(a,2U) ^ rotr32(a,13U) ^ rotr32(a,22U);
        uint32_t const maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t const t2 = s0 + maj;
        h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    ctx->state[0]+=a; ctx->state[1]+=b; ctx->state[2]+=c; ctx->state[3]+=d;
    ctx->state[4]+=e; ctx->state[5]+=f; ctx->state[6]+=g; ctx->state[7]+=h;
}

static void sha256_init(Sha256 *ctx)
{
    static const uint32_t init[8] = {
        0x6a09e667U,0xbb67ae85U,0x3c6ef372U,0xa54ff53aU,
        0x510e527fU,0x9b05688cU,0x1f83d9abU,0x5be0cd19U
    };
    memcpy(ctx->state, init, sizeof init);
    ctx->bitlen = 0U;
    ctx->datalen = 0U;
}

static void sha256_update(Sha256 *ctx, const void *ptr, size_t len)
{
    const uint8_t *data = ptr;
    for (size_t i = 0U; i < len; ++i) {
        ctx->data[ctx->datalen++] = data[i];
        if (ctx->datalen == 64U) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512U;
            ctx->datalen = 0U;
        }
    }
}

static void sha256_final(Sha256 *ctx, uint8_t out[32])
{
    size_t i = ctx->datalen;
    ctx->data[i++] = 0x80U;
    if (i > 56U) {
        while (i < 64U) ctx->data[i++] = 0U;
        sha256_transform(ctx, ctx->data);
        i = 0U;
    }
    while (i < 56U) ctx->data[i++] = 0U;
    ctx->bitlen += (uint64_t)ctx->datalen * 8U;
    for (unsigned j = 0U; j < 8U; ++j) {
        ctx->data[63U-j] = (uint8_t)(ctx->bitlen >> (8U*j));
    }
    sha256_transform(ctx, ctx->data);
    for (unsigned j = 0U; j < 8U; ++j) {
        out[4U*j]=(uint8_t)(ctx->state[j]>>24U);
        out[4U*j+1U]=(uint8_t)(ctx->state[j]>>16U);
        out[4U*j+2U]=(uint8_t)(ctx->state[j]>>8U);
        out[4U*j+3U]=(uint8_t)ctx->state[j];
    }
}

static void hash_to_hex(const uint8_t hash[32], char out[65])
{
    static const char hex[] = "0123456789abcdef";
    for (unsigned i = 0U; i < 32U; ++i) {
        out[2U*i] = hex[hash[i] >> 4U];
        out[2U*i+1U] = hex[hash[i] & 15U];
    }
    out[64] = '\0';
}

static bool sha256_file(const char *path, char out[65])
{
    int fd = open(path, O_RDONLY | O_CLOEXEC | O_NOFOLLOW);
    if (fd < 0) return false;
    Sha256 ctx;
    sha256_init(&ctx);
    uint8_t buf[65536];
    for (;;) {
        ssize_t const n = read(fd, buf, sizeof buf);
        if (n < 0) {
            if (errno == EINTR) continue;
            (void)close(fd);
            return false;
        }
        if (n == 0) break;
        sha256_update(&ctx, buf, (size_t)n);
    }
    if (close(fd) != 0) return false;
    uint8_t hash[32];
    sha256_final(&ctx, hash);
    hash_to_hex(hash, out);
    return true;
}

static bool sha256_current_executable(char out[65])
{
    char path[MAX_PATH_LEN];
    ssize_t const n=readlink("/proc/self/exe",path,sizeof path-1U);
    if(n<=0 || (size_t)n>=sizeof path) return false;
    path[n]='\0';
    return sha256_file(path,out);
}

static bool valid_hash(const char *s)
{
    if (strlen(s) != HASH_HEX_LEN) return false;
    for (size_t i = 0U; i < HASH_HEX_LEN; ++i) {
        if (!((s[i] >= '0' && s[i] <= '9') || (s[i] >= 'a' && s[i] <= 'f'))) return false;
    }
    return true;
}

static void path_join(char out[MAX_PATH_LEN], const char *dir, const char *base)
{
    int const n = snprintf(out, MAX_PATH_LEN, "%s/%s", dir, base);
    if (n < 0 || (unsigned)n >= MAX_PATH_LEN) die("path is too long");
}

static void write_all(int fd, const void *ptr, size_t len)
{
    const uint8_t *p = ptr;
    while (len != 0U) {
        ssize_t const n = write(fd, p, len);
        if (n < 0) {
            if (errno == EINTR) continue;
            die("write failed: %s", strerror(errno));
        }
        if (n == 0) die("short zero-byte write");
        p += (size_t)n;
        len -= (size_t)n;
    }
}

static void read_all(int fd, void *ptr, size_t len)
{
    uint8_t *p = ptr;
    while (len != 0U) {
        ssize_t const n = read(fd, p, len);
        if (n < 0) {
            if (errno == EINTR) continue;
            die("read failed: %s", strerror(errno));
        }
        if (n == 0) die("truncated file");
        p += (size_t)n;
        len -= (size_t)n;
    }
}

static void fsync_dir(const char *dir)
{
    int const fd = open(dir, O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW);
    if (fd < 0) die("cannot open run directory for fsync: %s", strerror(errno));
    if (fsync(fd) != 0) die("directory fsync failed: %s", strerror(errno));
    if (close(fd) != 0) die("directory close failed: %s", strerror(errno));
}

static void atomic_write(const char *dir, const char *base, const void *data, size_t len)
{
    char final_path[MAX_PATH_LEN], tmp_path[MAX_PATH_LEN];
    path_join(final_path, dir, base);
    int const n = snprintf(tmp_path, sizeof tmp_path, "%s/%s.tmp", dir, base);
    if (n < 0 || (size_t)n >= sizeof tmp_path) die("temporary path is too long");
    int const fd = open(tmp_path, O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC | O_NOFOLLOW, 0600);
    if (fd < 0) die("cannot create %s: %s", tmp_path, strerror(errno));
    write_all(fd, data, len);
    if (fsync(fd) != 0) die("fsync failed for %s: %s", tmp_path, strerror(errno));
    if (close(fd) != 0) die("close failed for %s: %s", tmp_path, strerror(errno));
    if (rename(tmp_path, final_path) != 0) die("rename to %s failed: %s", final_path, strerror(errno));
    fsync_dir(dir);
}

static void put_u32le(uint8_t *p, uint32_t x)
{
    p[0]=(uint8_t)x; p[1]=(uint8_t)(x>>8U); p[2]=(uint8_t)(x>>16U); p[3]=(uint8_t)(x>>24U);
}

static void put_u64le(uint8_t *p, uint64_t x)
{
    for (unsigned i=0U; i<8U; ++i) p[i]=(uint8_t)(x>>(8U*i));
}

static uint32_t get_u32le(const uint8_t *p)
{
    return (uint32_t)p[0] | ((uint32_t)p[1]<<8U) | ((uint32_t)p[2]<<16U) | ((uint32_t)p[3]<<24U);
}

static uint64_t get_u64le(const uint8_t *p)
{
    uint64_t x=0U;
    for (unsigned i=0U; i<8U; ++i) x |= (uint64_t)p[i] << (8U*i);
    return x;
}

static void build_binary_header(uint8_t h[BIN_HEADER_SIZE], uint32_t kind,
                                uint32_t B, uint32_t elem_size, uint64_t count,
                                const uint8_t payload_hash[32])
{
    memset(h, 0, BIN_HEADER_SIZE);
    memcpy(h, BINARY_MAGIC, strlen(BINARY_MAGIC));
    put_u32le(h+16U, FORMAT_VERSION);
    put_u32le(h+20U, kind);
    put_u32le(h+24U, B);
    put_u32le(h+28U, elem_size);
    put_u64le(h+32U, count);
    memcpy(h+40U, payload_hash, 24U); /* first 192 bits; whole-file SHA-256 is authoritative */
}

static void write_binary_atomic(const char *dir, const char *base, uint32_t kind,
                                uint32_t B, uint32_t elem_size, const void *payload,
                                uint64_t count)
{
    if (count > SIZE_MAX) die("binary element count exceeds size_t");
    size_t bytes;
    if (!mul_size((size_t)count, (size_t)elem_size, &bytes)) die("binary payload size overflow");
    Sha256 sh;
    uint8_t digest[32], header[BIN_HEADER_SIZE];
    sha256_init(&sh);
    sha256_update(&sh, payload, bytes);
    sha256_final(&sh, digest);
    build_binary_header(header, kind, B, elem_size, count, digest);

    char final_path[MAX_PATH_LEN], tmp_path[MAX_PATH_LEN];
    path_join(final_path, dir, base);
    int const n=snprintf(tmp_path,sizeof tmp_path,"%s/%s.tmp",dir,base);
    if (n<0 || (size_t)n>=sizeof tmp_path) die("temporary path is too long");
    int const fd=open(tmp_path,O_WRONLY|O_CREAT|O_EXCL|O_CLOEXEC|O_NOFOLLOW,0600);
    if (fd<0) die("cannot create %s: %s",tmp_path,strerror(errno));
    write_all(fd,header,sizeof header);
    write_all(fd,payload,bytes);
    if (fsync(fd)!=0) die("fsync failed for %s: %s",tmp_path,strerror(errno));
    if (close(fd)!=0) die("close failed for %s: %s",tmp_path,strerror(errno));
    test_fault_exit("before_rename",91);
    if (rename(tmp_path,final_path)!=0) die("rename to %s failed: %s",final_path,strerror(errno));
    test_fault_exit("after_rename",92);
    fsync_dir(dir);
}

static void *read_binary(const char *path, uint32_t expected_kind, uint32_t B,
                         uint32_t elem_size, uint64_t expected_count)
{
    int const fd=open(path,O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if (fd<0) die("cannot open %s: %s",path,strerror(errno));
    struct stat st;
    if (fstat(fd,&st)!=0 || !S_ISREG(st.st_mode)) die("invalid binary artifact %s",path);
    uint8_t h[BIN_HEADER_SIZE];
    read_all(fd,h,sizeof h);
    if (memcmp(h,BINARY_MAGIC,strlen(BINARY_MAGIC))!=0 ||
        get_u32le(h+16U)!=FORMAT_VERSION || get_u32le(h+20U)!=expected_kind ||
        get_u32le(h+24U)!=B || get_u32le(h+28U)!=elem_size ||
        get_u64le(h+32U)!=expected_count) die("binary header mismatch in %s",path);
    if (expected_count>SIZE_MAX) die("binary element count exceeds size_t");
    size_t bytes;
    if (!mul_size((size_t)expected_count,(size_t)elem_size,&bytes)) die("binary size overflow");
    if ((uint64_t)st.st_size != (uint64_t)BIN_HEADER_SIZE + (uint64_t)bytes) die("binary size mismatch in %s",path);
    void *payload=malloc(bytes==0U?1U:bytes);
    if (payload==NULL) die("allocation failed while reading %s",path);
    read_all(fd,payload,bytes);
    uint8_t extra;
    if (read(fd,&extra,1U)!=0) die("trailing data in %s",path);
    if (close(fd)!=0) die("close failed for %s",path);
    Sha256 sh; uint8_t digest[32];
    sha256_init(&sh); sha256_update(&sh,payload,bytes); sha256_final(&sh,digest);
    if (memcmp(h+40U,digest,24U)!=0) die("payload hash mismatch in %s",path);
    return payload;
}

static uint32_t parse_B(const char *s)
{
    if (s==NULL || *s=='\0') die("B is empty");
    for (const char *p=s; *p!='\0'; ++p) if (*p<'0'||*p>'9') die("B must be strict unsigned decimal");
    errno=0; char *end=NULL; unsigned long long const v=strtoull(s,&end,10);
    if (errno!=0 || end==s || *end!='\0' || v<2ULL || v>(unsigned long long)MAX_B)
        die("B must satisfy 2 <= B <= %u",MAX_B);
    return (uint32_t)v;
}

static double monotonic_seconds(void)
{
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC,&ts)!=0) die("clock_gettime failed");
    return (double)ts.tv_sec + (double)ts.tv_nsec/1.0e9;
}

static uint64_t monotonic_nanoseconds(void)
{
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC,&ts)!=0) die("clock_gettime failed");
    return (uint64_t)ts.tv_sec*1000000000ULL+(uint64_t)ts.tv_nsec;
}

static void heartbeat_maybe(Heartbeat *heartbeat)
{
    uint64_t const now=monotonic_nanoseconds();
    uint_fast64_t expected=atomic_load_explicit(&heartbeat->next_ns,memory_order_relaxed);
    if ((uint_fast64_t)now<expected) return;
    uint_fast64_t const replacement=(uint_fast64_t)now+HEARTBEAT_INTERVAL_NS;
    if (!atomic_compare_exchange_strong_explicit(&heartbeat->next_ns,&expected,replacement,
                                                  memory_order_relaxed,memory_order_relaxed))
        return;
    uint64_t const completed=(uint64_t)atomic_load_explicit(&heartbeat->completed,
                                                            memory_order_relaxed);
    double const elapsed=(double)(now-heartbeat->start_ns)/1.0e9;
    double const fraction=heartbeat->total==0U?1.0:(double)completed/(double)heartbeat->total;
    (void)fprintf(stderr,
        "heartbeat\tc_lo=%u\tc_hi=%u\tcompleted=%" PRIu64 "\ttotal=%" PRIu64
        "\tfraction=%.6f\telapsed=%.3f\n",
        heartbeat->c_lo,heartbeat->c_hi,completed,heartbeat->total,fraction,elapsed);
}

static uint32_t *build_spf(uint32_t B)
{
    size_t count=(size_t)B+1U, bytes;
    if (!mul_size(count,sizeof(uint32_t),&bytes)) die("SPF allocation size overflow");
    uint32_t *spf=calloc(1U,bytes);
    if (spf==NULL) die("cannot allocate SPF table");
    if (B>=1U) spf[1]=1U;
    for (uint32_t i=2U; i<=B; ++i) {
        if (spf[i]!=0U) continue;
        spf[i]=i;
        if ((uint64_t)i*(uint64_t)i>(uint64_t)B) continue;
        for (uint64_t j=(uint64_t)i*(uint64_t)i; j<=(uint64_t)B; j+=(uint64_t)i) {
            if (spf[(size_t)j]==0U) spf[(size_t)j]=i;
        }
    }
    return spf;
}

static unsigned factorize(uint32_t n, const uint32_t *spf, Factor f[MAX_FACTORS])
{
    unsigned nf=0U;
    while (n>1U) {
        if (nf>=MAX_FACTORS) die("factor-count bound exceeded");
        uint32_t const p=spf[n];
        if (p<2U || n%p!=0U) die("invalid SPF table at %u",n);
        uint32_t e=0U, q=1U;
        do {
            if (q>UINT32_MAX/p) die("prime-power overflow");
            q*=p; n/=p; ++e;
        } while (n>1U && spf[n]==p);
        f[nf].p=p; f[nf].e=e; f[nf].q=q;
        f[nf].h=0U; f[nf].modulus=0U; f[nf].generator=0U;
        ++nf;
    }
    return nf;
}

static unsigned distinct_prime_factors(uint32_t n, const uint32_t *spf,
                                       uint32_t out[MAX_FACTORS])
{
    unsigned k=0U;
    while (n>1U) {
        uint32_t const p=spf[n];
        if (p<2U || n%p!=0U || k>=MAX_FACTORS) die("invalid p-1 factorization");
        out[k++]=p;
        do n/=p; while (n>1U && spf[n]==p);
    }
    return k;
}

static uint64_t primitive_root_odd(uint32_t p, const uint32_t *spf, RootCache *cache)
{
    size_t const slot=(size_t)p & (ROOT_CACHE_SIZE-1U);
    if (cache->key[slot]==p) return cache->value[slot];
    uint32_t divisors[MAX_FACTORS];
    unsigned const nd=distinct_prime_factors(p-1U,spf,divisors);
    uint32_t g0=2U;
    for (;; ++g0) {
        if (g0>=p) die("primitive-root search failed for p=%u",p);
        bool ok=true;
        for (unsigned i=0U; i<nd; ++i) {
            if (pow_mod_u64(g0,(uint64_t)(p-1U)/divisors[i],p)==1U) { ok=false; break; }
        }
        if (ok) break;
    }
    uint64_t const p2=(uint64_t)p*(uint64_t)p;
    uint64_t g=(uint64_t)g0;
    if (pow_mod_u64(g,(uint64_t)p-1U,p2)==1U) g+=(uint64_t)p;
    if (pow_mod_u64(g,(uint64_t)p-1U,p2)==1U) die("primitive-root lift failed for p=%u",p);
    cache->key[slot]=p; cache->value[slot]=g;
    return g;
}

static void prepare_factor_metadata(uint32_t c, Factor f[MAX_FACTORS], unsigned nf)
{
    uint32_t const k=c-1U;
    for (unsigned i=0U; i<nf; ++i) {
        uint64_t const q=f[i].q;
        f[i].modulus=q*q;
        if (f[i].modulus==0U ||
            f[i].modulus>(uint64_t)MAX_B*(uint64_t)MAX_B)
            die("local modulus outside numeric ledger");
        if (f[i].p==2U) {
            f[i].h=1U;
            f[i].generator=1U;
        } else {
            f[i].h=gcd_u32(k,f[i].p-1U);
            f[i].generator=0U;
        }
    }
}

static void prepare_selected_generators(uint32_t c, Factor f[MAX_FACTORS], unsigned nf,
                                        unsigned selected, const uint32_t *spf,
                                        RootCache *cache)
{
    for (unsigned i=0U; i<nf; ++i) {
        if ((selected&(1U<<i))==0U || f[i].p==2U) continue;
        if (f[i].h==1U) {
            f[i].generator=1U;
            continue;
        }
        uint64_t const phi=f[i].modulus-f[i].modulus/f[i].p;
        uint64_t const primitive=primitive_root_odd(f[i].p,spf,cache);
        f[i].generator=pow_mod_u64(primitive,phi/f[i].h,f[i].modulus);
        if (pow_mod_u64(f[i].generator,f[i].h,f[i].modulus)!=1U)
            die("local generator order does not divide h for c=%u",c);
    }
}

static unsigned choose_factor_mask(uint32_t c, uint32_t B, const Factor *f, unsigned nf)
{
    uint64_t const c2=(uint64_t)c*(uint64_t)c;
    unsigned const full=(1U<<nf)-1U;
    if (c2<=(uint64_t)B) return full;
    uint64_t best_cost=UINT64_MAX;
    unsigned best=0U;
    for (unsigned mask=1U; mask<=full; ++mask) {
        uint64_t d=1U, cost=1U;
        for (unsigned i=0U; i<nf; ++i) if ((mask&(1U<<i))!=0U) {
            d*=f[i].q;
            if (cost>UINT32_MAX/f[i].h) { cost=UINT64_MAX; break; }
            cost*=f[i].h;
        }
        if (cost==UINT64_MAX || d*d<=(uint64_t)B) continue;
        if (cost<best_cost || (cost==best_cost && mask<best)) {
            best_cost=cost; best=mask;
        }
    }
    if (best==0U) die("no interval-aware factor subset for c=%u",c);
    return best;
}

static void metric_inc(uint64_t *x)
{
    if (*x==UINT64_MAX) die("telemetry counter overflow");
    ++*x;
}

static void count_inc(uint8_t *counts, uint32_t b, Metrics *m)
{
    metric_inc(&m->increments);
    if (counts[b]>=COUNTER_SAT-1U) {
        counts[b]=COUNTER_SAT;
    } else {
        ++counts[b];
    }
}

static bool omitted_factors_accept(uint64_t b, uint32_t c, const Factor *f,
                                   unsigned nf, unsigned selected, Metrics *m)
{
    uint64_t const k=(uint64_t)c-1U;
    for (unsigned i=0U; i<nf; ++i) {
        if ((selected&(1U<<i))!=0U) continue;
        if (f[i].p==2U) {
            metric_inc(&m->two_adic_tests);
            if (b%f[i].modulus!=1U) {
                metric_inc(&m->two_adic_rejects);
                return false;
            }
        } else {
            uint32_t const h=gcd_u32((uint32_t)k,f[i].p-1U);
            if (pow_mod_u64(b,h,f[i].p)!=1U) {
                metric_inc(&m->order_rejects);
                return false;
            }
            metric_inc(&m->valuation_tests);
            if (pow_mod_u64(b,(uint64_t)f[i].p-1U,f[i].modulus)!=1U) {
                metric_inc(&m->valuation_rejects);
                return false;
            }
        }
    }
    return true;
}

typedef struct {
    uint32_t c;
    uint32_t B;
    const Factor *f;
    unsigned nf;
    unsigned selected;
    uint8_t *counts;
    Metrics *metrics;
    Heartbeat *heartbeat;
} EnumContext;

static void accept_projected_root(const EnumContext *ctx, uint64_t r, uint64_t modulus)
{
    Metrics *m=ctx->metrics;
    metric_inc(&m->roots_generated);
    if ((m->roots_generated&(HEARTBEAT_ROOT_STRIDE-1U))==0U)
        heartbeat_maybe(ctx->heartbeat);
    uint64_t const c2=(uint64_t)ctx->c*(uint64_t)ctx->c;
    if (c2<=(uint64_t)ctx->B) {
        uint64_t b=r;
        if (b<=(uint64_t)ctx->c) {
            uint64_t const t=((uint64_t)ctx->c-b)/modulus+1U;
            b+=(uint64_t)((u128)t*(u128)modulus);
        }
        if (b>(uint64_t)ctx->B) {
            metric_inc(&m->range_rejects);
            return;
        }
        for (;;) {
            count_inc(ctx->counts,(uint32_t)b,m);
            if (modulus>(uint64_t)ctx->B-b) break;
            b+=modulus;
        }
    } else {
        if (r<=(uint64_t)ctx->c || r>(uint64_t)ctx->B) {
            metric_inc(&m->range_rejects);
            return;
        }
        if (omitted_factors_accept(r,ctx->c,ctx->f,ctx->nf,ctx->selected,m))
            count_inc(ctx->counts,(uint32_t)r,m);
    }
}

static void enumerate_roots_rec(const EnumContext *ctx, unsigned factor_index,
                                uint64_t x, uint64_t modulus)
{
    unsigned i=factor_index;
    while (i<ctx->nf && (ctx->selected&(1U<<i))==0U) ++i;
    if (i==ctx->nf) {
        accept_projected_root(ctx,x,modulus);
        return;
    }
    Factor const *fi=&ctx->f[i];
    uint64_t y=1U;
    uint64_t xm=0U, inv=0U;
    if (modulus!=1U) {
        xm=x%fi->modulus;
        inv=inv_mod_u64(modulus%fi->modulus,fi->modulus);
    }
    for (uint32_t j=0U; j<fi->h; ++j) {
        uint64_t nx, nm;
        if (modulus==1U) {
            nx=y; nm=fi->modulus;
        } else {
            uint64_t const diff=(y>=xm)?(y-xm):(fi->modulus-(xm-y));
            uint64_t const t=mul_mod_u64(diff,inv,fi->modulus);
            u128 const wide=(u128)x+(u128)modulus*(u128)t;
            u128 const wide_mod=(u128)modulus*(u128)fi->modulus;
            if (wide_mod>UINT64_MAX) die("CRT modulus overflow");
            nm=(uint64_t)wide_mod;
            nx=(uint64_t)(wide%wide_mod);
        }
        enumerate_roots_rec(ctx,i+1U,nx,nm);
        y=mul_mod_u64(y,fi->generator,fi->modulus);
    }
    if (y!=1U) die("local root enumeration did not close for c=%u",ctx->c);
}

static void process_composite(uint32_t c, uint32_t B, const uint32_t *spf,
                              RootCache *cache, uint8_t *counts, Metrics *metrics,
                              Heartbeat *heartbeat)
{
    Factor f[MAX_FACTORS];
    unsigned const nf=factorize(c,spf,f);
    if (nf==1U && f[0].e==1U) die("prime entered composite kernel");
    prepare_factor_metadata(c,f,nf);
    unsigned const selected=choose_factor_mask(c,B,f,nf);
    prepare_selected_generators(c,f,nf,selected,spf,cache);
    EnumContext const ctx={c,B,f,nf,selected,counts,metrics,heartbeat};
    metric_inc(&metrics->composites);
    enumerate_roots_rec(&ctx,0U,0U,1U);
}

static void metrics_add(Metrics *a, const Metrics *b)
{
#define ADD_FIELD(name) do { if (!add_u64(a->name,b->name,&a->name)) die("metric overflow: " #name); } while (0)
    ADD_FIELD(roots_generated); ADD_FIELD(range_rejects); ADD_FIELD(order_rejects);
    ADD_FIELD(valuation_tests); ADD_FIELD(valuation_rejects); ADD_FIELD(two_adic_tests);
    ADD_FIELD(two_adic_rejects); ADD_FIELD(increments); ADD_FIELD(saturations);
    ADD_FIELD(composites);
#undef ADD_FIELD
}

static unsigned choose_threads(uint32_t B)
{
    int requested=omp_get_max_threads();
    if (requested<1) requested=1;
    if (requested>64) requested=64;
    long const pages=sysconf(_SC_AVPHYS_PAGES), page_size=sysconf(_SC_PAGESIZE);
    uint64_t avail=0U;
    if (pages>0 && page_size>0) avail=(uint64_t)pages*(uint64_t)page_size;
    uint64_t const n=(uint64_t)B+1U;
    for (int t=requested; t>=1; --t) {
        /* SPF + active count + segment delta + T private counts + bounded scratch. */
        u128 const need=(u128)6U*n+(u128)(unsigned)t*n+(u128)64U*1024U*1024U;
        if (need>SIZE_MAX) continue;
        if (avail==0U || need*4U<=(u128)avail*3U) return (unsigned)t;
    }
    die("verified memory budget cannot support even one thread");
    return 1U;
}

static Metrics compute_segment(uint32_t B, uint32_t c_lo, uint32_t c_hi,
                               const uint32_t *spf, uint8_t *delta,
                               unsigned *threads_out)
{
    unsigned const requested_threads=choose_threads(B);
    size_t const count=(size_t)B+1U;
    if(requested_threads==0U || count<3U || count>SIZE_MAX-(CACHE_LINE_SIZE-1U))
        die("invalid thread/allocation domain");
    size_t const stride=(count+(CACHE_LINE_SIZE-1U))&~((size_t)CACHE_LINE_SIZE-1U);
    size_t elems, private_bytes, metrics_bytes, cache_bytes;
    if (!mul_size((size_t)requested_threads,stride,&elems) ||
        !mul_size(elems,sizeof(uint8_t),&private_bytes) ||
        !mul_size((size_t)requested_threads,sizeof(MetricsSlot),&metrics_bytes) ||
        !mul_size((size_t)requested_threads,sizeof(RootCacheSlot),&cache_bytes))
        die("thread-local allocation size overflow");
    uint8_t *private_counts=aligned_alloc(CACHE_LINE_SIZE,private_bytes);
    MetricsSlot *per=aligned_alloc(CACHE_LINE_SIZE,metrics_bytes);
    RootCacheSlot *caches=aligned_alloc(CACHE_LINE_SIZE,cache_bytes);
    if (private_counts==NULL || per==NULL || caches==NULL) die("thread-local allocation failed");
    memset(private_counts,0,private_bytes);
    memset(per,0,metrics_bytes);
    memset(caches,0,cache_bytes);

    Heartbeat heartbeat;
    uint64_t const heartbeat_start=monotonic_nanoseconds();
    atomic_init(&heartbeat.completed,0U);
    atomic_init(&heartbeat.next_ns,(uint_fast64_t)heartbeat_start+HEARTBEAT_INTERVAL_NS);
    heartbeat.start_ns=heartbeat_start;
    heartbeat.total=(uint64_t)c_hi-(uint64_t)c_lo+1U;
    heartbeat.c_lo=c_lo;
    heartbeat.c_hi=c_hi;
    unsigned actual_threads=0U;

#pragma omp parallel num_threads(requested_threads) default(none) shared(B,c_lo,c_hi,spf,private_counts,per,caches,stride,heartbeat,actual_threads)
    {
        unsigned const tid=(unsigned)omp_get_thread_num();
        uint8_t *local=private_counts+(size_t)tid*stride;
#pragma omp single
        actual_threads=(unsigned)omp_get_num_threads();
#pragma omp for schedule(dynamic,OMP_DYNAMIC_CHUNK)
        for (uint64_t cw=(uint64_t)c_lo; cw<=(uint64_t)c_hi; ++cw) {
            uint32_t const c=(uint32_t)cw;
            if (spf[c]!=c)
                process_composite(c,B,spf,&caches[tid].value,local,&per[tid].value,&heartbeat);
            uint64_t const completed=(uint64_t)atomic_fetch_add_explicit(
                &heartbeat.completed,1U,memory_order_relaxed)+1U;
            if (completed==heartbeat.total || completed%HEARTBEAT_COMPOSITE_STRIDE==0U)
                heartbeat_maybe(&heartbeat);
        }
    }

    if(actual_threads==0U || actual_threads>requested_threads)
        die("invalid actual OpenMP team size");
    *threads_out=actual_threads;

    Metrics total={0};
    memset(delta,0,count);
    for (unsigned t=0U; t<actual_threads; ++t) {
        metrics_add(&total,&per[t].value);
        uint8_t const *local=private_counts+(size_t)t*stride;
        for (uint32_t b=2U; b<=B; ++b) {
            unsigned const sum=(unsigned)delta[b]+(unsigned)local[b];
            delta[b]=(uint8_t)(sum>=COUNTER_SAT?COUNTER_SAT:sum);
        }
    }
    free(caches); free(per); free(private_counts);
    return total;
}

static FILE *open_text_read(const char *path)
{
    int const fd=open(path,O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if (fd<0) die("cannot open %s: %s",path,strerror(errno));
    struct stat st;
    if (fstat(fd,&st)!=0 || !S_ISREG(st.st_mode)) {
        (void)close(fd);
        die("invalid text artifact %s",path);
    }
    FILE *fp=fdopen(fd,"r");
    if (fp==NULL) { (void)close(fd); die("fdopen failed for %s",path); }
    return fp;
}

static uint32_t parse_u32_field(const char *s, const char *name)
{
    if (s==NULL || *s=='\0') die("empty %s",name);
    for (const char *p=s; *p!='\0'; ++p) if (*p<'0'||*p>'9') die("invalid %s",name);
    errno=0; char *end=NULL; unsigned long const v=strtoul(s,&end,10);
    if (errno!=0 || *end!='\0' || v>UINT32_MAX) die("out-of-range %s",name);
    return (uint32_t)v;
}

static uint64_t parse_u64_field(const char *s, const char *name)
{
    if (s==NULL || *s=='\0') die("empty %s",name);
    for (const char *p=s; *p!='\0'; ++p) if (*p<'0'||*p>'9') die("invalid %s",name);
    errno=0; char *end=NULL; unsigned long long const v=strtoull(s,&end,10);
    if (errno!=0 || *end!='\0') die("out-of-range %s",name);
    return (uint64_t)v;
}

static void write_config(const char *dir, const Config *cfg)
{
    char buf[2048];
    int const n=snprintf(buf,sizeof buf,
        "%s\t%u\n"
        "program_id\t%s\n"
        "B\t%u\n"
        "index_cap\t%u\n"
        "counter_saturation\t%u\n"
        "segment_span\t%u\n"
        "test_noncampaign\t%u\n"
        "source_sha256\t%s\n"
        "makefile_sha256\t%s\n"
        "binary_sha256\t%s\n"
        "artifact_count\t%u\n",
        CONFIG_MAGIC,FORMAT_VERSION,PROGRAM_ID,cfg->B,INDEX_CAP,COUNTER_SAT,
        cfg->segment_span,cfg->test_noncampaign,cfg->source_hash,cfg->makefile_hash,cfg->binary_hash,
        ARTIFACT_COUNT);
    if (n<0 || (size_t)n>=sizeof buf) die("config serialization overflow");
    atomic_write(dir,"config.tsv",buf,(size_t)n);
}

static Config read_config(const char *dir)
{
    char path[MAX_PATH_LEN]; path_join(path,dir,"config.tsv");
    FILE *fp=open_text_read(path);
    char *line=NULL; size_t cap=0U; ssize_t n;
    Config cfg={0}; unsigned seen=0U, lines=0U;
    while ((n=getline(&line,&cap,fp))>=0) {
        ++lines;
        if (n==0 || line[n-1]!='\n') die("unterminated config line");
        line[n-1]='\0';
        char *tab=strchr(line,'\t');
        if (tab==NULL || strchr(tab+1,'\t')!=NULL) die("malformed config line");
        *tab='\0'; const char *value=tab+1;
        if (lines==1U) {
            if (strcmp(line,CONFIG_MAGIC)!=0 || parse_u32_field(value,"config version")!=FORMAT_VERSION)
                die("config magic/version mismatch");
            continue;
        }
        unsigned bit=0U;
        if (strcmp(line,"program_id")==0) { bit=1U<<0U; if(strcmp(value,PROGRAM_ID)!=0) die("program_id mismatch"); }
        else if (strcmp(line,"B")==0) { bit=1U<<1U; cfg.B=parse_u32_field(value,"B"); }
        else if (strcmp(line,"index_cap")==0) { bit=1U<<2U; if(parse_u32_field(value,"index_cap")!=INDEX_CAP) die("index cap mismatch"); }
        else if (strcmp(line,"counter_saturation")==0) { bit=1U<<3U; if(parse_u32_field(value,"counter saturation")!=COUNTER_SAT) die("counter saturation mismatch"); }
        else if (strcmp(line,"segment_span")==0) { bit=1U<<4U; cfg.segment_span=parse_u32_field(value,"segment_span"); }
        else if (strcmp(line,"test_noncampaign")==0) { bit=1U<<5U; cfg.test_noncampaign=parse_u32_field(value,"test_noncampaign"); if(cfg.test_noncampaign>1U) die("invalid test_noncampaign"); }
        else if (strcmp(line,"source_sha256")==0) { bit=1U<<6U; if(!valid_hash(value)) die("invalid source hash"); memcpy(cfg.source_hash,value,65U); }
        else if (strcmp(line,"makefile_sha256")==0) { bit=1U<<7U; if(!valid_hash(value)) die("invalid Makefile hash"); memcpy(cfg.makefile_hash,value,65U); }
        else if (strcmp(line,"binary_sha256")==0) { bit=1U<<8U; if(!valid_hash(value)) die("invalid binary hash"); memcpy(cfg.binary_hash,value,65U); }
        else if (strcmp(line,"artifact_count")==0) { bit=1U<<9U; if(parse_u32_field(value,"artifact_count")!=ARTIFACT_COUNT) die("artifact count mismatch"); }
        else die("unknown config key: %s",line);
        if ((seen&bit)!=0U) die("duplicate config key: %s",line);
        seen|=bit;
    }
    free(line);
    if (ferror(fp)!=0 || fclose(fp)!=0) die("config read failure");
    if (lines!=11U || seen!=((1U<<10U)-1U)) die("incomplete config");
    uint32_t const work=(cfg.B>4U)?cfg.B-4U:0U;
    if (cfg.B<2U || cfg.B>MAX_B || cfg.segment_span==0U ||
        (work==0U?cfg.segment_span!=1U:cfg.segment_span>work))
        die("config numeric domain violation");
    if (cfg.test_noncampaign!=0U && cfg.B>2024U) die("non-campaign test bound violation");
    return cfg;
}

static void verify_build_identity(const Config *cfg)
{
    if (!valid_hash(SOURCE_SHA256) || !valid_hash(MAKEFILE_SHA256))
        die("binary was not built by the canonical Makefile");
    if (strcmp(cfg->source_hash,SOURCE_SHA256)!=0 || strcmp(cfg->makefile_hash,MAKEFILE_SHA256)!=0)
        die("source/Makefile identity mismatch");
    char current[65];
    if (!sha256_current_executable(current)) die("cannot hash current executable");
    if (strcmp(current,cfg->binary_hash)!=0) die("binary identity mismatch");
}

static const char *manifest_header(void)
{
    return "segment\tc_lo\tc_hi\told_frontier\tfrontier\troots_generated\trange_rejects\torder_rejects\tvaluation_tests\tvaluation_rejects\ttwo_adic_tests\ttwo_adic_rejects\tincrements\tsaturations\tcomposites\telapsed_seconds\tthreads\tstate\tactive_slot\tconfig_sha256\tspf_sha256\tcount_sha256\tsource_sha256\tbinary_sha256\tstop_reason\n";
}

static void serialize_manifest_row(const ManifestRow *r, const Config *cfg,
                                   char out[2048], size_t *len)
{
    int const n=snprintf(out,2048,
        "%" PRIu64 "\t%u\t%u\t%u\t%u\t"
        "%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t"
        "%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t"
        "%.9f\t%u\t%s\t%u\t%s\t%s\t%s\t%s\t%s\t%s\n",
        r->segment,r->c_lo,r->c_hi,r->old_frontier,r->frontier,
        r->metrics.roots_generated,r->metrics.range_rejects,r->metrics.order_rejects,
        r->metrics.valuation_tests,r->metrics.valuation_rejects,r->metrics.two_adic_tests,
        r->metrics.two_adic_rejects,r->metrics.increments,r->metrics.saturations,
        r->metrics.composites,r->elapsed,r->threads,r->state,r->active_slot,
        r->config_hash,r->spf_hash,r->count_hash,cfg->source_hash,cfg->binary_hash,
        r->stop_reason);
    if (n<0 || n>=2048) die("manifest row serialization overflow");
    *len=(size_t)n;
}

static void append_manifest(const char *dir, const ManifestRow *r, const Config *cfg)
{
    char path[MAX_PATH_LEN], row[2048]; size_t len;
    path_join(path,dir,"manifest.tsv");
    serialize_manifest_row(r,cfg,row,&len);
    int const fd=open(path,O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if(fd<0) die("cannot read manifest for append: %s",strerror(errno));
    struct stat st;
    if(fstat(fd,&st)!=0 || !S_ISREG(st.st_mode) || st.st_size<0 ||
       (uintmax_t)st.st_size>(uintmax_t)(SIZE_MAX-len)) {
        (void)close(fd);
        die("invalid manifest size before append");
    }
    size_t const old_len=(size_t)st.st_size;
    char *combined=malloc(old_len+len);
    if(combined==NULL) { (void)close(fd); die("manifest append allocation failed"); }
    read_all(fd,combined,old_len);
    if(close(fd)!=0) die("manifest close failed");
    memcpy(combined+old_len,row,len);
    atomic_write(dir,"manifest.tsv",combined,old_len+len);
    free(combined);
}

static unsigned split_tabs(char *line, char **field, unsigned max_fields)
{
    if(max_fields==0U) die("zero TSV field capacity");
    unsigned n=1U;
    field[0]=line;
    for(char *p=line; *p!='\0'; ++p) {
        if(*p!='\t') continue;
        *p='\0';
        if(n>=max_fields) die("too many TSV fields");
        field[n++]=p+1;
    }
    return n;
}

static bool metrics_are_zero(const Metrics *m)
{
    Metrics const zero={0};
    return memcmp(m,&zero,sizeof zero)==0;
}

static ManifestRow read_manifest(const char *dir, const Config *cfg, uint64_t *rows_out,
                                 ManifestRow *previous_out)
{
    char path[MAX_PATH_LEN]; path_join(path,dir,"manifest.tsv");
    char config_path[MAX_PATH_LEN], spf_path[MAX_PATH_LEN], actual_config[65], actual_spf[65];
    path_join(config_path,dir,"config.tsv"); path_join(spf_path,dir,"spf.bin");
    if(!sha256_file(config_path,actual_config) || !sha256_file(spf_path,actual_spf))
        die("cannot hash manifest-bound config/SPF");
    FILE *fp=open_text_read(path);
    char *line=NULL; size_t cap=0U; ssize_t nr;
    if ((nr=getline(&line,&cap,fp))<0 || nr==0 || line[nr-1]!='\n') die("truncated manifest magic");
    line[nr-1]='\0'; char expected_magic[64];
    (void)snprintf(expected_magic,sizeof expected_magic,"%s\t%u",MANIFEST_MAGIC,FORMAT_VERSION);
    if (strcmp(line,expected_magic)!=0) die("manifest magic mismatch");
    if (getline(&line,&cap,fp)<0 || strcmp(line,manifest_header())!=0)
        die("manifest header mismatch");
    ManifestRow last={0}, previous={0};
    uint64_t rows=0U; uint32_t expected_frontier=(cfg->B<4U)?cfg->B-1U:3U;
    while ((nr=getline(&line,&cap,fp))>=0) {
        if (nr==0 || line[nr-1]!='\n') die("truncated manifest row");
        line[nr-1]='\0';
        char *f[25]; unsigned const n=split_tabs(line,f,25U);
        if (n!=25U) die("manifest field-count mismatch");
        ManifestRow r={0};
        r.segment=parse_u64_field(f[0],"segment"); r.c_lo=parse_u32_field(f[1],"c_lo");
        r.c_hi=parse_u32_field(f[2],"c_hi"); r.old_frontier=parse_u32_field(f[3],"old_frontier");
        r.frontier=parse_u32_field(f[4],"frontier");
        r.metrics.roots_generated=parse_u64_field(f[5],"roots_generated");
        r.metrics.range_rejects=parse_u64_field(f[6],"range_rejects");
        r.metrics.order_rejects=parse_u64_field(f[7],"order_rejects");
        r.metrics.valuation_tests=parse_u64_field(f[8],"valuation_tests");
        r.metrics.valuation_rejects=parse_u64_field(f[9],"valuation_rejects");
        r.metrics.two_adic_tests=parse_u64_field(f[10],"two_adic_tests");
        r.metrics.two_adic_rejects=parse_u64_field(f[11],"two_adic_rejects");
        r.metrics.increments=parse_u64_field(f[12],"increments");
        r.metrics.saturations=parse_u64_field(f[13],"saturations");
        r.metrics.composites=parse_u64_field(f[14],"composites");
        errno=0; char *end=NULL; r.elapsed=strtod(f[15],&end);
        if(errno!=0 || end==f[15] || *end!='\0' || !isfinite(r.elapsed) || r.elapsed<0.0)
            die("invalid elapsed_seconds");
        r.threads=(unsigned)parse_u32_field(f[16],"threads");
        if(strlen(f[17])>=sizeof r.state) die("state too long");
        strcpy(r.state,f[17]);
        r.active_slot=(unsigned)parse_u32_field(f[18],"active_slot");
        for(unsigned hi=19U; hi<=23U; ++hi) if(!valid_hash(f[hi])) die("invalid manifest hash");
        strcpy(r.config_hash,f[19]); strcpy(r.spf_hash,f[20]); strcpy(r.count_hash,f[21]);
        if(strcmp(r.config_hash,actual_config)!=0 || strcmp(r.spf_hash,actual_spf)!=0 ||
           strcmp(f[22],cfg->source_hash)!=0 || strcmp(f[23],cfg->binary_hash)!=0)
            die("manifest identity mismatch");
        if(strlen(f[24])>=sizeof r.stop_reason) die("stop reason too long");
        strcpy(r.stop_reason,f[24]);
        if(r.segment!=rows || r.active_slot>1U || r.threads==0U || r.threads>64U)
            die("manifest sequence/slot/thread mismatch");
        if(r.old_frontier!=expected_frontier) die("noncontiguous manifest old frontier");
        if(rows==0U) {
            if(r.c_lo!=0U || r.c_hi!=0U || r.frontier!=expected_frontier ||
               r.active_slot!=0U || r.threads!=1U || r.elapsed!=0.0 ||
               !metrics_are_zero(&r.metrics) || strcmp(r.stop_reason,"initialized")!=0)
                die("invalid initial manifest row");
        } else {
            uint32_t const expected_lo=expected_frontier+1U;
            uint64_t expected_hi=(uint64_t)expected_lo+(uint64_t)cfg->segment_span-1U;
            if(expected_hi>(uint64_t)cfg->B-1U) expected_hi=(uint64_t)cfg->B-1U;
            if(r.c_lo!=expected_lo || r.c_hi!=(uint32_t)expected_hi ||
               r.frontier!=r.c_hi || r.frontier>=cfg->B ||
               r.active_slot!=1U-last.active_slot)
                die("noncontiguous manifest coverage");
            bool const complete=r.frontier>=cfg->B-1U;
            if(strcmp(r.stop_reason,complete?"complete":"segment_limit")!=0)
                die("manifest stop-reason/frontier mismatch");
            previous=last;
        }
        bool const complete=r.frontier>=cfg->B-1U;
        if(strcmp(r.state,complete?"complete":"incomplete")!=0) die("manifest state/frontier mismatch");
        expected_frontier=r.frontier;
        last=r;
        ++rows;
    }
    free(line);
    if (ferror(fp)!=0 || fclose(fp)!=0 || rows==0U) die("manifest read failure");
    *rows_out=rows;
    if(previous_out!=NULL) *previous_out=previous;
    return last;
}

static size_t build_results_text(char buf[RESULTS_BUFFER_SIZE], const Config *cfg,
                                 const ManifestRow *last, const uint8_t *counts)
{
    size_t used=0U;
    bool const complete=last->frontier>=cfg->B-1U;
    uint32_t base_frontier=cfg->B;
    uint32_t minima[INDEX_CAP+1U];
    memset(minima,0,sizeof minima);
    if(complete) {
        for(uint32_t b=2U; b<=cfg->B; ++b) {
            uint8_t const target=counts[b];
            if(target>=1U && target<=INDEX_CAP && minima[target]==0U)
                minima[target]=b;
        }
    }
    if (!complete && cfg->B>=4U) base_frontier=last->frontier+1U;
    int n=snprintf(buf+used,RESULTS_BUFFER_SIZE-used,
        "A255885-RESULTS\t%u\nstate\t%s\nB\t%u\nmodulus_frontier\t%u\nbase_frontier\t%u\n"
        "index_cap\t%u\ncounter_saturation\t%u\nindex\tbase\n",
        FORMAT_VERSION,complete?"complete":"incomplete",cfg->B,last->frontier,
        base_frontier,INDEX_CAP,COUNTER_SAT);
    if(n<0 || (size_t)n>=RESULTS_BUFFER_SIZE-used) die("results serialization overflow");
    used+=(size_t)n;
    for(uint32_t target=1U; target<=INDEX_CAP; ++target) {
        uint32_t const found=minima[target];
        n=snprintf(buf+used,RESULTS_BUFFER_SIZE-used,"%u\t%s",target,found==0U?"NA":"");
        if(n<0 || (size_t)n>=RESULTS_BUFFER_SIZE-used) die("results serialization overflow");
        used+=(size_t)n;
        if(found!=0U) {
            n=snprintf(buf+used,RESULTS_BUFFER_SIZE-used,"%u",found);
            if(n<0 || (size_t)n>=RESULTS_BUFFER_SIZE-used) die("results serialization overflow");
            used+=(size_t)n;
        }
        if(used+1U>=RESULTS_BUFFER_SIZE) die("results serialization overflow");
        buf[used++]='\n';
    }
    return used;
}

static void write_results(const char *dir, const Config *cfg, const ManifestRow *last,
                          const uint8_t *counts)
{
    char buf[RESULTS_BUFFER_SIZE];
    size_t const used=build_results_text(buf,cfg,last,counts);
    atomic_write(dir,"results.tsv",buf,used);
}

static const char *const checksum_files[6]={
    "config.tsv","cnt0.bin","cnt1.bin","manifest.tsv","results.tsv","spf.bin"
};

typedef struct {
    char hash[6][HASH_HEX_LEN+1U];
} ChecksumSnapshot;

static void write_checksums(const char *dir)
{
    char buf[2048]; size_t used=0U;
    for(unsigned i=0U; i<6U; ++i) {
        char path[MAX_PATH_LEN], hash[65]; path_join(path,dir,checksum_files[i]);
        if(!sha256_file(path,hash)) die("cannot hash %s",path);
        int const n=snprintf(buf+used,sizeof buf-used,"%s  %s\n",hash,checksum_files[i]);
        if(n<0 || (size_t)n>=sizeof buf-used) die("checksum serialization overflow");
        used+=(size_t)n;
    }
    atomic_write(dir,"checksums.sha256",buf,used);
}

static bool read_checksum_snapshot(const char *dir, ChecksumSnapshot *snapshot)
{
    char path[MAX_PATH_LEN]; path_join(path,dir,"checksums.sha256");
    int const fd=open(path,O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if(fd<0) return false;
    struct stat st;
    if(fstat(fd,&st)!=0 || !S_ISREG(st.st_mode)) { (void)close(fd); return false; }
    FILE *fp=fdopen(fd,"r");
    if(fp==NULL) { (void)close(fd); return false; }
    char *line=NULL; size_t cap=0U; ssize_t n;
    unsigned row=0U;
    while((n=getline(&line,&cap,fp))>=0) {
        if(row>=6U || n!=(ssize_t)(64U+2U+strlen(checksum_files[row])+1U) || line[n-1]!='\n')
            { free(line); (void)fclose(fp); return false; }
        line[n-1]='\0';
        if(line[64]!=' ' || line[65]!=' ' || strcmp(line+66,checksum_files[row])!=0)
            { free(line); (void)fclose(fp); return false; }
        line[64]='\0';
        if(!valid_hash(line)) { free(line); (void)fclose(fp); return false; }
        memcpy(snapshot->hash[row],line,HASH_HEX_LEN+1U);
        ++row;
    }
    free(line);
    bool const ok=ferror(fp)==0 && row==6U;
    if(fclose(fp)!=0) return false;
    return ok;
}

static bool file_hash_matches(const char *dir, const char *base, const char *expected)
{
    char path[MAX_PATH_LEN], actual[HASH_HEX_LEN+1U];
    path_join(path,dir,base);
    return sha256_file(path,actual) && strcmp(actual,expected)==0;
}

static bool checksums_match(const char *dir)
{
    ChecksumSnapshot snapshot;
    if(!read_checksum_snapshot(dir,&snapshot)) return false;
    for(unsigned i=0U; i<6U; ++i)
        if(!file_hash_matches(dir,checksum_files[i],snapshot.hash[i])) return false;
    return true;
}

static bool manifest_prefix_hash(const char *dir, char out[HASH_HEX_LEN+1U])
{
    char path[MAX_PATH_LEN]; path_join(path,dir,"manifest.tsv");
    int const fd=open(path,O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if(fd<0) return false;
    struct stat st;
    if(fstat(fd,&st)!=0 || !S_ISREG(st.st_mode) || st.st_size<=0 ||
       (uintmax_t)st.st_size>SIZE_MAX) { (void)close(fd); return false; }
    size_t const len=(size_t)st.st_size;
    uint8_t *data=malloc(len);
    if(data==NULL) { (void)close(fd); return false; }
    read_all(fd,data,len);
    if(close(fd)!=0) { free(data); return false; }
    if(data[len-1U]!='\n') { free(data); return false; }
    size_t cut=len-1U;
    while(cut>0U && data[cut-1U]!='\n') --cut;
    if(cut==0U) { free(data); return false; }
    Sha256 sh; uint8_t digest[32];
    sha256_init(&sh); sha256_update(&sh,data,cut); sha256_final(&sh,digest);
    hash_to_hex(digest,out);
    free(data);
    return true;
}

static void verify_checksums(const char *dir)
{
    if(!checksums_match(dir)) die("checksum verification failed");
}

static bool allowed_artifact(const char *name)
{
    static const char *const names[ARTIFACT_COUNT]={
        "config.tsv","spf.bin","cnt0.bin","cnt1.bin","manifest.tsv","results.tsv","checksums.sha256"
    };
    for(unsigned i=0U;i<ARTIFACT_COUNT;++i) if(strcmp(name,names[i])==0) return true;
    return false;
}

static void verify_artifact_set(const char *dir)
{
    DIR *dp=opendir(dir); if(dp==NULL) die("cannot open run directory: %s",strerror(errno));
    unsigned count=0U; struct dirent *de;
    errno=0;
    while((de=readdir(dp))!=NULL) {
        if(strcmp(de->d_name,".")==0 || strcmp(de->d_name,"..")==0) continue;
        if(!allowed_artifact(de->d_name)) die("unexpected run artifact: %s",de->d_name);
        ++count;
    }
    if(errno!=0 || closedir(dp)!=0) die("run directory read failure");
    if(count!=ARTIFACT_COUNT) die("run artifact count mismatch: %u",count);
}

static int lock_run_directory(const char *dir, int operation)
{
    int const fd=open(dir,O_RDONLY|O_DIRECTORY|O_CLOEXEC|O_NOFOLLOW);
    if(fd<0) die("cannot open run directory for locking: %s",strerror(errno));
    while(flock(fd,operation)!=0) {
        int const saved_errno=errno;
        if(saved_errno==EINTR) continue;
        (void)close(fd);
        if(saved_errno==EWOULDBLOCK || saved_errno==EAGAIN)
            die("run directory is already in use by another a255885 process");
        die("cannot lock run directory: %s",strerror(saved_errno));
    }
    return fd;
}

static void release_run_directory_lock(int fd)
{
    (void)close(fd);
}

static bool test_stop_after_init_lock_requested(const Config *cfg)
{
    char const *stop=getenv("A255885_TEST_STOP_AFTER_INIT_LOCK");
    if(stop==NULL) return false;
    if(strcmp(stop,"1")!=0 || cfg->test_noncampaign==0U)
        die("init-lock stop is restricted to recorded non-campaign tests");
    return true;
}

static void test_stop_after_init_lock(bool requested)
{
    if(requested && raise(SIGSTOP)!=0)
        die("cannot enter init-lock test stop: %s",strerror(errno));
}

static void discard_recovery_temps(const char *dir)
{
    static const char *const names[]={
        "config.tsv.tmp","spf.bin.tmp","cnt0.bin.tmp","cnt1.bin.tmp",
        "manifest.tsv.tmp","results.tsv.tmp","checksums.sha256.tmp"
    };
    bool removed=false;
    for(size_t i=0U; i<sizeof names/sizeof names[0]; ++i) {
        char path[MAX_PATH_LEN]; path_join(path,dir,names[i]);
        if(unlink(path)==0) {
            removed=true;
            (void)fprintf(stderr,"recovery\tdiscarded=%s\n",names[i]);
        } else if(errno!=ENOENT) {
            die("cannot discard recovery temp %s: %s",names[i],strerror(errno));
        }
    }
    if(removed) fsync_dir(dir);
}

static void verify_results(const char *dir, const Config *cfg, const ManifestRow *last,
                           const uint8_t *counts)
{
    char expected[RESULTS_BUFFER_SIZE], path[MAX_PATH_LEN];
    size_t const expected_len=build_results_text(expected,cfg,last,counts);
    path_join(path,dir,"results.tsv");
    int const fd=open(path,O_RDONLY|O_CLOEXEC|O_NOFOLLOW);
    if(fd<0) die("cannot open results.tsv");
    struct stat st;
    if(fstat(fd,&st)!=0 || !S_ISREG(st.st_mode) || (uint64_t)st.st_size!=(uint64_t)expected_len)
        die("results.tsv size mismatch");
    char *actual=malloc(expected_len==0U?1U:expected_len);
    if(actual==NULL) die("results verification allocation failed");
    read_all(fd,actual,expected_len);
    if(close(fd)!=0) die("results close failed");
    if(memcmp(actual,expected,expected_len)!=0) die("results.tsv semantic mismatch");
    free(actual);
}

static void verify_spf_table(uint32_t B, const uint32_t *spf)
{
    uint32_t *expected=build_spf(B);
    size_t bytes;
    if(!mul_size((size_t)B+1U,sizeof(uint32_t),&bytes))
        die("SPF verification size overflow");
    if(memcmp(spf,expected,bytes)!=0)
        die("SPF table differs from deterministic reconstruction");
    free(expected);
}

static uint8_t *load_counts(const char *dir, unsigned slot, uint32_t B)
{
    char path[MAX_PATH_LEN];
    path_join(path,dir,slot==0U?"cnt0.bin":"cnt1.bin");
    return read_binary(path,2U,B,1U,(uint64_t)B+1U);
}

static uint32_t *load_spf(const char *dir, uint32_t B)
{
    char path[MAX_PATH_LEN]; path_join(path,dir,"spf.bin");
    return read_binary(path,1U,B,4U,(uint64_t)B+1U);
}

static void verify_active_hash(const char *dir, const ManifestRow *last)
{
    char path[MAX_PATH_LEN], actual[65];
    path_join(path,dir,last->active_slot==0U?"cnt0.bin":"cnt1.bin");
    if(!sha256_file(path,actual) || strcmp(actual,last->count_hash)!=0)
        die("active checkpoint does not match manifest");
}

static void recover_interrupted_commit(const char *dir, const Config *cfg,
                                       const ManifestRow *last,
                                       const ManifestRow *previous, uint64_t rows)
{
    ChecksumSnapshot snapshot;
    if(!read_checksum_snapshot(dir,&snapshot))
        die("checksum ledger is malformed; recovery refused");
    if(!file_hash_matches(dir,"config.tsv",snapshot.hash[0]) ||
       !file_hash_matches(dir,"spf.bin",snapshot.hash[5]))
        die("unrecognized immutable-artifact mismatch; recovery refused");

    bool const current_manifest=file_hash_matches(dir,"manifest.tsv",snapshot.hash[3]);
    char prefix_hash[HASH_HEX_LEN+1U];
    bool const appended_manifest=!current_manifest && rows>=2U &&
        manifest_prefix_hash(dir,prefix_hash) && strcmp(prefix_hash,snapshot.hash[3])==0;

    if(current_manifest) {
        unsigned const active=last->active_slot;
        unsigned const inactive=1U-active;
        if(!file_hash_matches(dir,"results.tsv",snapshot.hash[4]) ||
           !file_hash_matches(dir,active==0U?"cnt0.bin":"cnt1.bin",snapshot.hash[1U+active]) ||
           file_hash_matches(dir,inactive==0U?"cnt0.bin":"cnt1.bin",snapshot.hash[1U+inactive]))
            die("checksum mismatch is not an interrupted checkpoint rename");
        uint8_t *authoritative=load_counts(dir,active,cfg->B);
        write_binary_atomic(dir,inactive==0U?"cnt0.bin":"cnt1.bin",2U,cfg->B,1U,
                            authoritative,(uint64_t)cfg->B+1U);
        write_checksums(dir);
        free(authoritative);
        (void)fprintf(stderr,
            "recovery\tstate=discarded_uncommitted_checkpoint\tsegment=%" PRIu64
            "\tactive_slot=%u\n",last->segment,last->active_slot);
        return;
    }

    if(appended_manifest) {
        unsigned const old_active=previous->active_slot;
        if(!file_hash_matches(dir,old_active==0U?"cnt0.bin":"cnt1.bin",
                              snapshot.hash[1U+old_active]) ||
           !file_hash_matches(dir,old_active==0U?"cnt0.bin":"cnt1.bin",
                              previous->count_hash))
            die("manifest prefix does not bind the previous checkpoint");
        verify_active_hash(dir,last);
        uint8_t *authoritative=load_counts(dir,last->active_slot,cfg->B);
        write_results(dir,cfg,last,authoritative);
        write_checksums(dir);
        free(authoritative);
        (void)fprintf(stderr,
            "recovery\tstate=completed_manifest_commit\tsegment=%" PRIu64
            "\tactive_slot=%u\n",last->segment,last->active_slot);
        return;
    }

    die("checksum mismatch has no recognized interrupted-commit lineage");
}

static void initialize_run(const char *dir, uint32_t B, const char *segment_arg)
{
    if(!valid_hash(SOURCE_SHA256) || !valid_hash(MAKEFILE_SHA256))
        die("use the canonical Makefile to build this program");
    Config cfg={0};
    cfg.B=B;
    uint32_t const work=(B>4U)?B-4U:0U;
    cfg.segment_span=work==0U?1U:(work+15U)/16U;
    char const *test_span=getenv("A255885_TEST_SEGMENT_SPAN");
    char const *noncampaign=getenv("A255885_NONCAMPAIGN");
    if(segment_arg!=NULL && (test_span!=NULL || noncampaign!=NULL))
        die("explicit segment span cannot be combined with test overrides");
    if(segment_arg!=NULL) {
        uint32_t const span=parse_u32_field(segment_arg,"SEGMENT_SPAN");
        if(span==0U || (work==0U?span!=1U:span>work))
            die("segment span outside work range");
        cfg.segment_span=span;
    } else if(test_span!=NULL) {
        if(noncampaign==NULL || strcmp(noncampaign,"1")!=0 || B>2024U)
            die("test segment override requires A255885_NONCAMPAIGN=1 and B<=2024");
        uint32_t const span=parse_u32_field(test_span,"A255885_TEST_SEGMENT_SPAN");
        if(span==0U || (work!=0U && span>work)) die("test segment span outside work range");
        cfg.segment_span=span;
        cfg.test_noncampaign=1U;
    } else if(noncampaign!=NULL) {
        die("A255885_NONCAMPAIGN is accepted only with A255885_TEST_SEGMENT_SPAN");
    }
    bool const stop_after_lock=test_stop_after_init_lock_requested(&cfg);
    if(mkdir(dir,0700)!=0) die("init refuses existing/uncreatable RUN_DIR: %s",strerror(errno));
    int const lock_fd=lock_run_directory(dir,LOCK_EX);
    test_stop_after_init_lock(stop_after_lock);
    fsync_dir(dir);
    strcpy(cfg.source_hash,SOURCE_SHA256);
    strcpy(cfg.makefile_hash,MAKEFILE_SHA256);
    if(!sha256_current_executable(cfg.binary_hash)) die("cannot hash current executable");
    write_config(dir,&cfg);

    uint32_t *spf=build_spf(B);
    size_t const count=(size_t)B+1U;
    uint8_t *counts=calloc(count,1U);
    if(counts==NULL) die("count allocation failed");
    write_binary_atomic(dir,"spf.bin",1U,B,4U,spf,(uint64_t)count);
    write_binary_atomic(dir,"cnt0.bin",2U,B,1U,counts,(uint64_t)count);
    write_binary_atomic(dir,"cnt1.bin",2U,B,1U,counts,(uint64_t)count);

    char config_path[MAX_PATH_LEN], spf_path[MAX_PATH_LEN], count_path[MAX_PATH_LEN];
    char config_hash[65], spf_hash[65], count_hash[65];
    path_join(config_path,dir,"config.tsv"); path_join(spf_path,dir,"spf.bin");
    path_join(count_path,dir,"cnt0.bin");
    if(!sha256_file(config_path,config_hash) || !sha256_file(spf_path,spf_hash) ||
       !sha256_file(count_path,count_hash)) die("cannot hash initial state");
    ManifestRow row={0};
    row.segment=0U; row.c_lo=0U; row.c_hi=0U;
    row.frontier=(B<4U)?B-1U:3U; row.old_frontier=row.frontier;
    row.threads=1U; row.active_slot=0U; row.elapsed=0.0;
    strcpy(row.state,row.frontier>=B-1U?"complete":"incomplete");
    strcpy(row.config_hash,config_hash); strcpy(row.spf_hash,spf_hash);
    strcpy(row.count_hash,count_hash); strcpy(row.stop_reason,"initialized");

    char manifest[4096], row_text[2048]; size_t row_len;
    serialize_manifest_row(&row,&cfg,row_text,&row_len);
    int const n=snprintf(manifest,sizeof manifest,"%s\t%u\n%s",MANIFEST_MAGIC,FORMAT_VERSION,manifest_header());
    if(n<0 || (size_t)n+row_len>=sizeof manifest) die("initial manifest serialization overflow");
    memcpy(manifest+(size_t)n,row_text,row_len);
    atomic_write(dir,"manifest.tsv",manifest,(size_t)n+row_len);
    write_results(dir,&cfg,&row,counts);
    write_checksums(dir);
    verify_artifact_set(dir);
    free(counts); free(spf);
    (void)fprintf(stdout,"initialized\t%s\tB=%u\tsegment_span=%u\tstate=%s\n",
                  dir,B,cfg.segment_span,row.state);
    release_run_directory_lock(lock_fd);
}

static void status_run(const char *dir)
{
    int const lock_fd=lock_run_directory(dir,LOCK_SH|LOCK_NB);
    verify_artifact_set(dir);
    Config const cfg=read_config(dir);
    verify_build_identity(&cfg);
    verify_checksums(dir);
    uint64_t rows=0U; ManifestRow const last=read_manifest(dir,&cfg,&rows,NULL);
    verify_active_hash(dir,&last);
    uint32_t const base_frontier=(last.frontier>=cfg.B-1U || cfg.B<4U)?cfg.B:last.frontier+1U;
    (void)fprintf(stdout,
        "program_id\t%s\nB\t%u\nstate\t%s\nsegments_committed\t%" PRIu64
        "\nmodulus_frontier\t%u\nbase_frontier\t%u\nactive_slot\t%u\nstop_reason\t%s\n",
        PROGRAM_ID,cfg.B,last.state,rows-1U,last.frontier,base_frontier,
        last.active_slot,last.stop_reason);
    release_run_directory_lock(lock_fd);
}

static void next_segment(const char *dir)
{
    int const lock_fd=lock_run_directory(dir,LOCK_EX|LOCK_NB);
    discard_recovery_temps(dir);
    verify_artifact_set(dir);
    Config const cfg=read_config(dir);
    verify_build_identity(&cfg);
    ManifestRow previous={0};
    uint64_t rows=0U; ManifestRow const last=read_manifest(dir,&cfg,&rows,&previous);
    verify_active_hash(dir,&last);
    const char *fault=getenv("A255885_TEST_FAULT");
    if(fault!=NULL) {
        bool const known=strcmp(fault,"before_rename")==0 ||
                         strcmp(fault,"after_rename")==0 ||
                         strcmp(fault,"after_manifest")==0 ||
                         strcmp(fault,"after_results")==0;
        if(!known || cfg.test_noncampaign==0U)
            die("fault injection is restricted to recorded non-campaign tests");
    }
    if(!checksums_match(dir)) {
        recover_interrupted_commit(dir,&cfg,&last,&previous,rows);
    }
    verify_checksums(dir);
    if(last.frontier>=cfg.B-1U) {
        (void)fprintf(stdout,"already_complete\tB=%u\n",cfg.B);
        release_run_directory_lock(lock_fd);
        return;
    }
    uint32_t const c_lo=last.frontier+1U;
    uint64_t hi64=(uint64_t)c_lo+(uint64_t)cfg.segment_span-1U;
    if(hi64>(uint64_t)cfg.B-1U) hi64=(uint64_t)cfg.B-1U;
    uint32_t const c_hi=(uint32_t)hi64;
    uint32_t *spf=load_spf(dir,cfg.B);
    uint8_t *counts=load_counts(dir,last.active_slot,cfg.B);
    uint8_t *delta=calloc((size_t)cfg.B+1U,1U);
    if(delta==NULL) die("segment delta allocation failed");

    (void)fprintf(stderr,
        "start\tB=%u\tindex_cap=%u\tsegment=%" PRIu64 "\tc_lo=%u\tc_hi=%u\tfrontier=%u\n",
        cfg.B,INDEX_CAP,rows,c_lo,c_hi,last.frontier);
    double const start=monotonic_seconds();
    unsigned threads=1U;
    Metrics metrics=compute_segment(cfg.B,c_lo,c_hi,spf,delta,&threads);
    for(uint32_t b=2U; b<=cfg.B; ++b) {
        unsigned const sum=(unsigned)counts[b]+(unsigned)delta[b];
        if(sum>=COUNTER_SAT) {
            if(counts[b]<COUNTER_SAT) metric_inc(&metrics.saturations);
            counts[b]=COUNTER_SAT;
        } else counts[b]=(uint8_t)sum;
    }
    double const elapsed=monotonic_seconds()-start;
    unsigned const inactive=1U-last.active_slot;
    g_test_fault_armed=cfg.test_noncampaign!=0U;
    write_binary_atomic(dir,inactive==0U?"cnt0.bin":"cnt1.bin",2U,cfg.B,1U,
                        counts,(uint64_t)cfg.B+1U);
    g_test_fault_armed=false;
    char cp_path[MAX_PATH_LEN], cp_hash[65];
    path_join(cp_path,dir,inactive==0U?"cnt0.bin":"cnt1.bin");
    if(!sha256_file(cp_path,cp_hash)) die("cannot hash new checkpoint");

    ManifestRow row={0};
    row.segment=rows; row.c_lo=c_lo; row.c_hi=c_hi;
    row.old_frontier=last.frontier; row.frontier=c_hi; row.threads=threads;
    row.active_slot=inactive; row.metrics=metrics; row.elapsed=elapsed;
    strcpy(row.state,c_hi>=cfg.B-1U?"complete":"incomplete");
    strcpy(row.config_hash,last.config_hash); strcpy(row.spf_hash,last.spf_hash);
    strcpy(row.count_hash,cp_hash);
    strcpy(row.stop_reason,c_hi>=cfg.B-1U?"complete":"segment_limit");
    append_manifest(dir,&row,&cfg);
    g_test_fault_armed=cfg.test_noncampaign!=0U;
    test_fault_exit("after_manifest",93);
    g_test_fault_armed=false;
    write_results(dir,&cfg,&row,counts);
    g_test_fault_armed=cfg.test_noncampaign!=0U;
    test_fault_exit("after_results",94);
    g_test_fault_armed=false;
    write_checksums(dir);
    (void)fprintf(stderr,
        "close\tB=%u\tsegment=%" PRIu64 "\tfrontier=%u\troots=%" PRIu64
        "\trange_rejects=%" PRIu64 "\torder_rejects=%" PRIu64
        "\tvaluation_tests=%" PRIu64 "\tvaluation_rejects=%" PRIu64
        "\ttwo_adic_tests=%" PRIu64 "\ttwo_adic_rejects=%" PRIu64
        "\tincrements=%" PRIu64 "\tsaturations=%" PRIu64
        "\telapsed=%.6f\tthreads=%u\tstop_reason=%s\n",
        cfg.B,row.segment,row.frontier,metrics.roots_generated,metrics.range_rejects,
        metrics.order_rejects,metrics.valuation_tests,metrics.valuation_rejects,
        metrics.two_adic_tests,metrics.two_adic_rejects,metrics.increments,
        metrics.saturations,elapsed,threads,row.stop_reason);
    free(delta); free(counts); free(spf);
    release_run_directory_lock(lock_fd);
}

static void verify_run(const char *dir)
{
    int const lock_fd=lock_run_directory(dir,LOCK_SH|LOCK_NB);
    verify_artifact_set(dir);
    Config const cfg=read_config(dir);
    verify_build_identity(&cfg);
    verify_checksums(dir);
    uint64_t rows=0U; ManifestRow const last=read_manifest(dir,&cfg,&rows,NULL);
    verify_active_hash(dir,&last);
    uint32_t *spf=load_spf(dir,cfg.B);
    uint8_t *active=load_counts(dir,last.active_slot,cfg.B);
    uint8_t *inactive=load_counts(dir,1U-last.active_slot,cfg.B);
    verify_spf_table(cfg.B,spf);
#if COUNTER_SAT < UINT8_MAX
    for(uint32_t b=0U; b<=cfg.B; ++b) {
        if(active[b]>COUNTER_SAT || inactive[b]>COUNTER_SAT) die("counter domain violation at %u",b);
    }
#endif
    verify_results(dir,&cfg,&last,active);
    uint32_t const base_frontier=(last.frontier>=cfg.B-1U || cfg.B<4U)?cfg.B:last.frontier+1U;
    (void)fprintf(stdout,
        "VERIFY PASS\tscope=identity+persistence+coverage\tarithmetic=not_recomputed"
        "\tstate=%s\tB=%u\tmodulus_frontier=%u\tbase_frontier=%u"
        "\tmanifest_rows=%" PRIu64 "\tactive_slot=%u\n",
        last.state,cfg.B,last.frontier,base_frontier,rows,last.active_slot);
    free(inactive); free(active); free(spf);
    release_run_directory_lock(lock_fd);
}

static void usage(FILE *stream)
{
    (void)fprintf(stream,
        "Usage:\n"
        "  a255885 init RUN_DIR B [SEGMENT_SPAN]\n"
        "  a255885 status RUN_DIR\n"
        "  a255885 next RUN_DIR\n"
        "  a255885 verify RUN_DIR\n"
        "\nB is strict decimal in [2,%u]. SEGMENT_SPAN is optional and"
        " covers composite-modulus work units per next invocation.\n"
        "init never overwrites RUN_DIR.\n",MAX_B);
}

int main(int argc, char **argv)
{
    if(argc<2) { usage(stderr); return EXIT_FAILURE; }
    if(strcmp(argv[1],"init")==0) {
        if(argc!=4 && argc!=5) { usage(stderr); return EXIT_FAILURE; }
        initialize_run(argv[2],parse_B(argv[3]),argc==5?argv[4]:NULL);
    } else if(strcmp(argv[1],"status")==0) {
        if(argc!=3) { usage(stderr); return EXIT_FAILURE; }
        status_run(argv[2]);
    } else if(strcmp(argv[1],"next")==0) {
        if(argc!=3) { usage(stderr); return EXIT_FAILURE; }
        next_segment(argv[2]);
    } else if(strcmp(argv[1],"verify")==0) {
        if(argc!=3) { usage(stderr); return EXIT_FAILURE; }
        verify_run(argv[2]);
    } else {
        usage(stderr);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
