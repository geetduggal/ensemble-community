// Simple fixed-memory linear probing hash table
// Geet Duggal

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
  // Basic properties
  uint32_t size_proper;
  uint16_t size_key;
  uint16_t size_entry;

  // Handy registers
  uint8_t log2size;
  uint8_t * last_entry;
  uint32_t modmask;
 
  uint8_t * entries;
  uint32_t num_collisions;
  uint32_t num_entries;
} ht_t;

// Memory management
ht_t * ht_alloc(uint8_t log2size, uint16_t size_key, uint16_t size_entry);
void ht_free(ht_t * ht);

// Basic API
uint8_t * ht_find(ht_t * ht, uint8_t * key);
uint8_t * ht_insert(ht_t * ht, uint8_t * found, uint8_t * entry);
void ht_delete(ht_t * ht, uint8_t * found);

// Test print functions
static inline uint8_t ht_occupied(ht_t * ht, uint8_t * entry);
void ht_printblock(ht_t * ht);
uint32_t ht_contigram(ht_t * ht);

// Read/write to disk (simple for now)
uint32_t ht_fwrite(ht_t * ht, FILE * f);
ht_t * ht_fread(FILE * f);

// Is an entry occupied? 
static inline uint8_t ht_occupied(ht_t * ht, uint8_t * entry) {
  uint16_t size;
  
  size = ht->size_entry;
  do {
    size--;
    if (entry[size] != 255)
      return 1;
  } while (size > 0);

  return 0;
}

// FNV hash function
static inline uint32_t hash (uint8_t * key, uint16_t size_key) {
  uint32_t fnv;
  uint8_t * ekey;
 
  ekey = key + size_key;
  fnv = 0;
  while (key < ekey) {
    fnv ^= (uint32_t)*key++; 
    fnv += (fnv<<1) + (fnv<<4) + (fnv<<7) + (fnv<<8) + (fnv<<24);
  }

  return fnv;
}
