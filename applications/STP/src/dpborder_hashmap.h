/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dpborder_hashmap.h
 * @brief  Special hashmap for DP-border Steiner tree algorithm
 * @author Daniel Rehfeldt
 *
 * modified version of https://github.com/sheredom/hashmap.h
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_DPBORDER_HASHMAP_H_
#define APPLICATIONS_STP_SRC_DPBORDER_HASHMAP_H_

#if defined(_MSC_VER)
// Workaround a bug in the MSVC runtime where it uses __cplusplus when not
// defined.
#pragma warning(push, 0)
#pragma warning(disable : 4668)
#endif
#include <stdlib.h>
#include <string.h>
#include "scip/scip.h"


#if (defined(_MSC_VER) && defined(__AVX__)) ||                                 \
    (!defined(_MSC_VER) && defined(__SSE4_2__)) ||                            \
    (defined(__INTEL_COMPILER))
#define HASHMAP_SSE42
#endif

#if defined(HASHMAP_SSE42)
#include <nmmintrin.h>
#endif

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

#if defined(_MSC_VER)
#pragma warning(push)
/* Stop MSVC complaining about unreferenced functions */
#pragma warning(disable : 4505)
/* Stop MSVC complaining about not inlining functions */
#pragma warning(disable : 4710)
/* Stop MSVC complaining about inlining functions! */
#pragma warning(disable : 4711)
#elif defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif

#if defined(_MSC_VER)
#define HASHMAP_USED
#elif defined(__GNUC__)
#define HASHMAP_USED __attribute__((used))
#else
#define HASHMAP_USED
#endif

/* We need to keep keys and values. */
struct hashmap_element_s {
  unsigned key_len;
  int key_position;
  int in_use;
  int value;
};

/* A hashmap has some maximum size and current size, as well as the data to
 * hold. */
typedef struct hashmap_s {
  const char* keyarray; /* non owned! */
  unsigned table_size;
  unsigned size;
  struct hashmap_element_s *elements;
} DPBHASHMAP;

#define HASHMAP_MAX_CHAIN_LENGTH (8)

#if defined(__cplusplus)
extern "C" {
#endif

/// @brief Create a hashmap.
/// @param initial_size The initial size of the hashmap. Must be a power of two.
/// @param out_hashmap The storage for the created hashmap.
/// @param keyarr Array with keys
/// @return On success 0 is returned.
///
/// Note that the initial size of the hashmap must be a power of two, and
/// creation of the hashmap will fail if this is not the case.
static SCIP_RETCODE hashmap_create(const unsigned initial_size,
                          const char* keyarr,
                          struct hashmap_s *const out_hashmap) HASHMAP_USED;

static void hashmap_updateKeyarr(
                   const char* keyarr,
                   struct hashmap_s* out_hashmap) HASHMAP_USED;

/// @brief Put an element into the hashmap.
/// @param hashmap The hashmap to insert into.
/// @param key_position The string key position to use.
/// @param len The length of the string key.
/// @param value The value to insert.
/// @return On success 0 is returned.
///
/// The key string slice is not copied when creating the hashmap entry, and thus
/// must remain a valid pointer until the hashmap entry is removed or the
/// hashmap is destroyed.
static int hashmap_putcheck(struct hashmap_s *const hashmap, int key_position,
         const unsigned len, int value) HASHMAP_USED;

static void hashmap_put(struct hashmap_s *const hashmap, int position,
                       const unsigned len, int value) HASHMAP_USED;

/// @brief Get an element from the hashmap.
/// @param hashmap The hashmap to get from.
/// @param position The string key position to use.
/// @param len The length of the string key.
/// @return The previously set element, or NULL if none exists.
static int hashmap_get(const struct hashmap_s *const hashmap, int position,
                         const unsigned len) HASHMAP_USED;

/// @brief Remove an element from the hashmap.
/// @param hashmap The hashmap to remove from.
/// @param position The string key position to use.
/// @param len The length of the string key.
/// @return On success 0 is returned.
static int hashmap_remove(struct hashmap_s *const hashmap, int position,
                  const unsigned len) HASHMAP_USED;

/// @brief Iterate over all the elements in a hashmap.
/// @param hashmap The hashmap to iterate over.
/// @param hashmap_new The new hashmap.
/// @return If the entire hashmap was iterated then 0 is returned.
static int hashmap_iterate_pairs(struct hashmap_s *const hashmap,
                                 struct hashmap_s *const hashmap_new) HASHMAP_USED;

/// @brief Get the size of the hashmap.
/// @param hashmap The hashmap to get the size of.
/// @return The size of the hashmap.
static unsigned
hashmap_num_entries(const struct hashmap_s *const hashmap) HASHMAP_USED;

static SCIP_Bool
hashmap_isEmpty(const struct hashmap_s *const hashmap) HASHMAP_USED;

/// @brief Destroy the hashmap.
/// @param hashmap The hashmap to destroy.
static void hashmap_destroy(struct hashmap_s *const hashmap) HASHMAP_USED;

static unsigned hashmap_crc32_helper(const char *const s,
                                     const unsigned len) HASHMAP_USED;
static unsigned hashmap_hash_helper_int_helper(const struct hashmap_s *const m,
                                               const char *const keystring,
                                               const unsigned len) HASHMAP_USED;
static int hashmap_hash_helper(const struct hashmap_s *const m,
                               int position, const unsigned len,
                               unsigned *const out_index) HASHMAP_USED;
static int
hashmap_rehash_iterator(struct hashmap_s *const m_new,
                        struct hashmap_element_s *const e) HASHMAP_USED;
static int hashmap_rehash_helper(struct hashmap_s *const m) HASHMAP_USED;

#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)
#define HASHMAP_CAST(type, x) static_cast<type>(x)
#define HASHMAP_PTR_CAST(type, x) reinterpret_cast<type>(x)
#define HASHMAP_NULL NULL
#else
#define HASHMAP_CAST(type, x) ((type)x)
#define HASHMAP_PTR_CAST(type, x) ((type)x)
#define HASHMAP_NULL 0
#endif

static inline
int hashmap_match_helper(const struct hashmap_element_s *const element,
                         const char *const key, const char *const keyarr, const unsigned len) {
  return (element->key_len == len) && (0 == memcmp(&keyarr[element->key_position], key, len));
}


SCIP_RETCODE hashmap_create(const unsigned initial_size,
                   const char* keyarr,
                   struct hashmap_s *const out_hashmap) {
  out_hashmap->table_size = initial_size;
  out_hashmap->size = 0;
  out_hashmap->keyarray = keyarr;

  if (0 == initial_size || 0 != (initial_size & (initial_size - 1))) {
    return SCIP_ERROR;
  }

  out_hashmap->elements =
      HASHMAP_CAST(struct hashmap_element_s *,
                   calloc(initial_size, sizeof(struct hashmap_element_s)));
  if (!out_hashmap->elements) {
    return SCIP_ERROR;
  }

  return SCIP_OKAY;
}


void hashmap_updateKeyarr(
                   const char* keyarr,
                   struct hashmap_s* out_hashmap)
{
   assert(keyarr && out_hashmap);
   out_hashmap->keyarray = keyarr;
}


int hashmap_putcheck(struct hashmap_s *const m, int key_position,
                const unsigned len, int value) {
  unsigned int index;
  assert(key_position >= 0);

  /* Find a place to put our value. */
  while (!hashmap_hash_helper(m, key_position, len, &index)) {
    if (hashmap_rehash_helper(m)) {
      return 1;
    }
  }

  /* Set the data. */
  m->elements[index].key_position = key_position;
  m->elements[index].key_len = len;
  m->elements[index].value = value;

  /* If the hashmap element was not already in use, set that it is being used
   * and bump our size. */
  if (0 == m->elements[index].in_use) {
    m->elements[index].in_use = 1;
    m->size++;
  }

  return 0;
}


void hashmap_put(struct hashmap_s *const m, int position,
               const unsigned len, int value) {
  unsigned int index;
  assert(position >= 0);

  /* Find a place to put our value. */
  while (!hashmap_hash_helper(m, position, len, &index)) {
    if (hashmap_rehash_helper(m)) {
      assert(0);
    }
  }

  /* Set the data. */
  m->elements[index].key_position = position;
  m->elements[index].key_len = len;
  m->elements[index].value = value;

  /* If the hashmap element was not already in use, set that it is being used
   * and bump our size. */
  if (0 == m->elements[index].in_use) {
    m->elements[index].in_use = 1;
    m->size++;
  }
}

int hashmap_get(const struct hashmap_s *const m, int position,
         const unsigned len) {
  unsigned int curr;
  unsigned int i;
  const char* const key = &(m->keyarray[position]);
  const struct hashmap_element_s* const elements = m->elements;

  assert(m->keyarray);

  /* Find data location */
  curr = hashmap_hash_helper_int_helper(m, key, len);

  /* Linear probing, if necessary */
  for (i = 0; i < HASHMAP_MAX_CHAIN_LENGTH; i++) {
    if (elements[curr].in_use) {
      if (hashmap_match_helper(&elements[curr], key, m->keyarray, len)) {
        return elements[curr].value;
      }
    }

    curr = (curr + 1) % m->table_size;
  }

  /* Not found */
  return -1;
}

int hashmap_remove(struct hashmap_s *const m, int position,
      const unsigned len) {
  unsigned int i;
  unsigned int curr;
  const char* const key = &(m->keyarray[position]);

  /* Find key */
  curr = hashmap_hash_helper_int_helper(m, key, len);

  /* Linear probing, if necessary */
  for (i = 0; i < HASHMAP_MAX_CHAIN_LENGTH; i++) {
    if (m->elements[curr].in_use) {
      if (hashmap_match_helper(&m->elements[curr], key, m->keyarray, len)) {
        /* Blank out the fields including in_use */
        memset(&m->elements[curr], 0, sizeof(struct hashmap_element_s));

        /* Reduce the size */
        m->size--;

        return 0;
      }
    }

    curr = (curr + 1) % m->table_size;
  }

  return 1;
}

int hashmap_iterate_pairs(struct hashmap_s *const hashmap,
                          struct hashmap_s *const hashmap_new) {
  unsigned int i;
  struct hashmap_element_s *p;
  int r;

  /* Linear probing */
  for (i = 0; i < hashmap->table_size; i++) {
    p = &hashmap->elements[i];
    if (p->in_use) {
      r = hashmap_rehash_iterator(hashmap_new, p);
      switch (r) {
      case -1: /* remove item */
        memset(p, 0, sizeof(struct hashmap_element_s));
        hashmap->size--;
        break;
      case 0: /* continue iterating */
        break;
      default: /* early exit */
        return 1;
      }
    }
  }
  return 0;
}

void hashmap_destroy(struct hashmap_s *const m) {
  free(m->elements);
  memset(m, 0, sizeof(struct hashmap_s));
}

unsigned hashmap_num_entries(const struct hashmap_s *const m) {
  return m->size;
}


SCIP_Bool hashmap_isEmpty(const struct hashmap_s *const m) {
  return (m->size == 0);
}

unsigned hashmap_crc32_helper(const char *const s, const unsigned len) {
  unsigned i;
  unsigned crc32val = 0;

#if defined(HASHMAP_SSE42)
  for (i = 0; i < len; i++) {
    crc32val = _mm_crc32_u8(crc32val, HASHMAP_CAST(unsigned char, s[i]));
  }

  return crc32val;
#else
  // Using polynomial 0x11EDC6F41 to match SSE 4.2's crc function.
  static const unsigned crc32_tab[] = {
      0x00000000U, 0xF26B8303U, 0xE13B70F7U, 0x1350F3F4U, 0xC79A971FU,
      0x35F1141CU, 0x26A1E7E8U, 0xD4CA64EBU, 0x8AD958CFU, 0x78B2DBCCU,
      0x6BE22838U, 0x9989AB3BU, 0x4D43CFD0U, 0xBF284CD3U, 0xAC78BF27U,
      0x5E133C24U, 0x105EC76FU, 0xE235446CU, 0xF165B798U, 0x030E349BU,
      0xD7C45070U, 0x25AFD373U, 0x36FF2087U, 0xC494A384U, 0x9A879FA0U,
      0x68EC1CA3U, 0x7BBCEF57U, 0x89D76C54U, 0x5D1D08BFU, 0xAF768BBCU,
      0xBC267848U, 0x4E4DFB4BU, 0x20BD8EDEU, 0xD2D60DDDU, 0xC186FE29U,
      0x33ED7D2AU, 0xE72719C1U, 0x154C9AC2U, 0x061C6936U, 0xF477EA35U,
      0xAA64D611U, 0x580F5512U, 0x4B5FA6E6U, 0xB93425E5U, 0x6DFE410EU,
      0x9F95C20DU, 0x8CC531F9U, 0x7EAEB2FAU, 0x30E349B1U, 0xC288CAB2U,
      0xD1D83946U, 0x23B3BA45U, 0xF779DEAEU, 0x05125DADU, 0x1642AE59U,
      0xE4292D5AU, 0xBA3A117EU, 0x4851927DU, 0x5B016189U, 0xA96AE28AU,
      0x7DA08661U, 0x8FCB0562U, 0x9C9BF696U, 0x6EF07595U, 0x417B1DBCU,
      0xB3109EBFU, 0xA0406D4BU, 0x522BEE48U, 0x86E18AA3U, 0x748A09A0U,
      0x67DAFA54U, 0x95B17957U, 0xCBA24573U, 0x39C9C670U, 0x2A993584U,
      0xD8F2B687U, 0x0C38D26CU, 0xFE53516FU, 0xED03A29BU, 0x1F682198U,
      0x5125DAD3U, 0xA34E59D0U, 0xB01EAA24U, 0x42752927U, 0x96BF4DCCU,
      0x64D4CECFU, 0x77843D3BU, 0x85EFBE38U, 0xDBFC821CU, 0x2997011FU,
      0x3AC7F2EBU, 0xC8AC71E8U, 0x1C661503U, 0xEE0D9600U, 0xFD5D65F4U,
      0x0F36E6F7U, 0x61C69362U, 0x93AD1061U, 0x80FDE395U, 0x72966096U,
      0xA65C047DU, 0x5437877EU, 0x4767748AU, 0xB50CF789U, 0xEB1FCBADU,
      0x197448AEU, 0x0A24BB5AU, 0xF84F3859U, 0x2C855CB2U, 0xDEEEDFB1U,
      0xCDBE2C45U, 0x3FD5AF46U, 0x7198540DU, 0x83F3D70EU, 0x90A324FAU,
      0x62C8A7F9U, 0xB602C312U, 0x44694011U, 0x5739B3E5U, 0xA55230E6U,
      0xFB410CC2U, 0x092A8FC1U, 0x1A7A7C35U, 0xE811FF36U, 0x3CDB9BDDU,
      0xCEB018DEU, 0xDDE0EB2AU, 0x2F8B6829U, 0x82F63B78U, 0x709DB87BU,
      0x63CD4B8FU, 0x91A6C88CU, 0x456CAC67U, 0xB7072F64U, 0xA457DC90U,
      0x563C5F93U, 0x082F63B7U, 0xFA44E0B4U, 0xE9141340U, 0x1B7F9043U,
      0xCFB5F4A8U, 0x3DDE77ABU, 0x2E8E845FU, 0xDCE5075CU, 0x92A8FC17U,
      0x60C37F14U, 0x73938CE0U, 0x81F80FE3U, 0x55326B08U, 0xA759E80BU,
      0xB4091BFFU, 0x466298FCU, 0x1871A4D8U, 0xEA1A27DBU, 0xF94AD42FU,
      0x0B21572CU, 0xDFEB33C7U, 0x2D80B0C4U, 0x3ED04330U, 0xCCBBC033U,
      0xA24BB5A6U, 0x502036A5U, 0x4370C551U, 0xB11B4652U, 0x65D122B9U,
      0x97BAA1BAU, 0x84EA524EU, 0x7681D14DU, 0x2892ED69U, 0xDAF96E6AU,
      0xC9A99D9EU, 0x3BC21E9DU, 0xEF087A76U, 0x1D63F975U, 0x0E330A81U,
      0xFC588982U, 0xB21572C9U, 0x407EF1CAU, 0x532E023EU, 0xA145813DU,
      0x758FE5D6U, 0x87E466D5U, 0x94B49521U, 0x66DF1622U, 0x38CC2A06U,
      0xCAA7A905U, 0xD9F75AF1U, 0x2B9CD9F2U, 0xFF56BD19U, 0x0D3D3E1AU,
      0x1E6DCDEEU, 0xEC064EEDU, 0xC38D26C4U, 0x31E6A5C7U, 0x22B65633U,
      0xD0DDD530U, 0x0417B1DBU, 0xF67C32D8U, 0xE52CC12CU, 0x1747422FU,
      0x49547E0BU, 0xBB3FFD08U, 0xA86F0EFCU, 0x5A048DFFU, 0x8ECEE914U,
      0x7CA56A17U, 0x6FF599E3U, 0x9D9E1AE0U, 0xD3D3E1ABU, 0x21B862A8U,
      0x32E8915CU, 0xC083125FU, 0x144976B4U, 0xE622F5B7U, 0xF5720643U,
      0x07198540U, 0x590AB964U, 0xAB613A67U, 0xB831C993U, 0x4A5A4A90U,
      0x9E902E7BU, 0x6CFBAD78U, 0x7FAB5E8CU, 0x8DC0DD8FU, 0xE330A81AU,
      0x115B2B19U, 0x020BD8EDU, 0xF0605BEEU, 0x24AA3F05U, 0xD6C1BC06U,
      0xC5914FF2U, 0x37FACCF1U, 0x69E9F0D5U, 0x9B8273D6U, 0x88D28022U,
      0x7AB90321U, 0xAE7367CAU, 0x5C18E4C9U, 0x4F48173DU, 0xBD23943EU,
      0xF36E6F75U, 0x0105EC76U, 0x12551F82U, 0xE03E9C81U, 0x34F4F86AU,
      0xC69F7B69U, 0xD5CF889DU, 0x27A40B9EU, 0x79B737BAU, 0x8BDCB4B9U,
      0x988C474DU, 0x6AE7C44EU, 0xBE2DA0A5U, 0x4C4623A6U, 0x5F16D052U,
      0xAD7D5351U};

  for (i = 0; i < len; i++) {
    crc32val = crc32_tab[(HASHMAP_CAST(unsigned char, crc32val) ^
                          HASHMAP_CAST(unsigned char, s[i]))] ^
               (crc32val >> 8);
  }
  return crc32val;
#endif
}

unsigned hashmap_hash_helper_int_helper(const struct hashmap_s *const m,
                                        const char *const keystring,
                                        const unsigned len) {
  unsigned key = hashmap_crc32_helper(keystring, len);

  /* Robert Jenkins' 32 bit Mix Function */
  key += (key << 12);
  key ^= (key >> 22);
  key += (key << 4);
  key ^= (key >> 9);
  key += (key << 10);
  key ^= (key >> 2);
  key += (key << 7);
  key ^= (key >> 12);

  /* Knuth's Multiplicative Method */
  key = (key >> 3) * 2654435761;

  return key % m->table_size;
}

int hashmap_hash_helper(const struct hashmap_s *const m, int position,
                        const unsigned len, unsigned *const out_index) {
  unsigned int start, curr;
  unsigned int i;
  int total_in_use;
  const char* const key = &(m->keyarray[position]);

  /* If full, return immediately */
  if (m->size >= m->table_size) {
    return 0;
  }

  /* Find the best index */
  curr = start = hashmap_hash_helper_int_helper(m, key, len);

  /* First linear probe to check if we've already insert the element */
  total_in_use = 0;

  for (i = 0; i < HASHMAP_MAX_CHAIN_LENGTH; i++) {
    const int in_use = m->elements[curr].in_use;

    total_in_use += in_use;

    if (in_use && hashmap_match_helper(&m->elements[curr], key, m->keyarray, len)) {
      *out_index = curr;
      return 1;
    }

    curr = (curr + 1) % m->table_size;
  }

  curr = start;

  /* Second linear probe to actually insert our element (only if there was at
   * least one empty entry) */
  if (HASHMAP_MAX_CHAIN_LENGTH > total_in_use) {
    for (i = 0; i < HASHMAP_MAX_CHAIN_LENGTH; i++) {
      if (!m->elements[curr].in_use) {
        *out_index = curr;
        return 1;
      }

      curr = (curr + 1) % m->table_size;
    }
  }

  return 0;
}

int hashmap_rehash_iterator(struct hashmap_s *const new_hash,
                            struct hashmap_element_s *const e) {
  const int temp = hashmap_putcheck(new_hash, e->key_position,
                         e->key_len, e->value);
  if (0 < temp) {
    return 1;
  }
  /* clear old value to avoid stale pointers */
  return -1;
}

/*
 * Doubles the size of the hashmap, and rehashes all the elements
 */
int hashmap_rehash_helper(struct hashmap_s *const m) {
  /* If this multiplication overflows hashmap_create will fail. */
  unsigned new_size = 2 * m->table_size;
  struct hashmap_s new_hash;
  SCIP_RETCODE retcode = hashmap_create(new_size, m->keyarray, &new_hash);
  int flag;

  printf("reallocating hash map \n");

  if (SCIP_OKAY != retcode) {
    return 1;
  }

  /* copy the old elements to the new table */
  flag = hashmap_iterate_pairs(m, &new_hash);

  if (0 != flag) {
    return flag;
  }

  hashmap_destroy(m);
  /* put new hash into old hash structure by copying */
  memcpy(m, &new_hash, sizeof(struct hashmap_s));

  return 0;
}

#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__clang__)
#pragma clang diagnostic pop
#endif



#endif /* APPLICATIONS_STP_SRC_DPBORDER_HASHMAP_H_ */
