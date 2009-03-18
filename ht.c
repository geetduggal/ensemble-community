#include "ht.h"
#include <string.h>

ht_t * ht_alloc(uint8_t log2size, uint16_t size_key, uint16_t size_entry) {
  uint32_t i;
  if (size_key > size_entry || size_entry < 1) 
    return NULL;

  ht_t * ht = malloc(sizeof(ht_t));
  
  ht->log2size = log2size;
  ht->size_proper = (1<<log2size);
  ht->size_key = size_key;
  ht->size_entry = size_entry;
  
  ht->num_entries = 0;
  
  ht->entries = calloc(ht->size_proper, size_entry);
  for (i = 0; i<ht->size_proper*size_entry; i++)
    ht->entries[i] = 255;

  ht->modmask = ht->size_proper-1;
  ht->last_entry = ht->entries + ht->size_entry*ht->modmask;

  return ht;
}

void ht_free(ht_t * ht) {
  free(ht->entries);
  free(ht);;
}

uint8_t * ht_find(ht_t * ht, uint8_t * key) {
  uint32_t num_probes;
  uint8_t * found;

  // The natural hash index
  found = &(ht->entries[(hash(key, ht->size_key)&ht->modmask)*ht->size_entry]);

  // Loop through the circular array
  num_probes = 0;
  do {
    // If we've found or entry or an unoccupied entry
    if ((memcmp(found, key, ht->size_key) == 0)
      || !ht_occupied(ht, found))
      return found;

    // Circular array without mod
    if (found == ht->last_entry)
      found = ht->entries;
    else
      found += ht->size_entry;

    num_probes++;
  } while (num_probes < ht->size_proper);

  // Looped through all items in the table
  return NULL;
} 

uint8_t * ht_insert(ht_t * ht, uint8_t * found, uint8_t * entry) {
  // We are full
  if (ht->num_entries == ht->size_proper)
    return NULL;
  // Simply copy the contents of the entry into the location
  else
    memcpy(found, entry, ht->size_entry);

  ht->num_entries++;
  return found;
}

void ht_delete(ht_t * ht, uint8_t * found) {
  uint32_t num_probes;
  uint8_t * next, * next_natural;

  // Get the found hash index from the pointer  
  next = found;
  
  num_probes = 0;

  do {
    // Circular array without mod
   if (next == ht->last_entry)
     next = ht->entries;
   else
     next += ht->size_entry;

   // Case 1: no other entry to worry about
   if (!ht_occupied(ht, next))
     break;
  
   // Case 2: uh oh, we have an entry to worry about
   next_natural = 
        &ht->entries[(hash(next, ht->size_key)&ht->modmask)*ht->size_entry];

   if ((next > found && (next_natural <= found || next_natural > next))
    || (next < found && (next_natural <= found && next_natural > next)) ) {
     // Move next entry back to fill in space
     memcpy(found, next, ht->size_entry);
     found = next;
   }

   num_probes++;
  } while (num_probes < ht->size_proper);

  // Clear proper entry
  memset(found, '\0', ht->size_entry);
}

void ht_printblock(ht_t * ht) {
  uint8_t * e;
  uint32_t i;

  for (e = ht->entries; e <= ht->last_entry; e += ht->size_entry) {
    i = (e - ht->entries)/ht->size_entry + 1;
    
    if (ht_occupied(ht, e))
      printf("*");
    else
      printf("-");

    if (i%78 == 0)
      printf("\n");
  }
}

uint32_t ht_contigram(ht_t * ht) {
  uint8_t o1, o2;
  uint32_t cnt, num_regions, num_clusters, size_cluster;
  uint8_t * entry;

  entry = ht->entries;
  // If the first slot is filled
  if (ht_occupied(ht, entry))
    printf("Starts filled\n");
  else
    printf("Starts empty\n");

  cnt = 1;
  num_regions = size_cluster = num_clusters = 0;
  // For all slots in the table
  while (entry <= ht->last_entry) {
    // If i and i+1 are either filled or empty
    o1 = ht_occupied(ht, ht->entries); 
    o2 = ht_occupied(ht, ht->entries+ht->size_entry); 

    if ((o1 && o2) || (!o1 && !o2))
      cnt++;
    else {
      if (o1 && !o2) {
        num_regions++;
        printf("%3u\t", cnt);
        if (num_regions%10 == 0)
          printf("\n");

        if (cnt > 3) {
          size_cluster += cnt;
          num_clusters++;
        }
      }
      cnt = 1;
    }
    entry += ht->size_entry;
  }
  return size_cluster/num_clusters;
}

uint32_t ht_fwrite(ht_t * ht, FILE * f) {
  uint32_t num_bytes;
  num_bytes = fwrite(ht, sizeof(ht_t), 1, f);
  num_bytes += fwrite(ht->entries, ht->size_entry, ht->size_proper, f);
  return num_bytes;
}

ht_t * ht_fread(FILE * f) {
  ht_t meta, * ht;
  fread(&meta, sizeof(ht_t), 1, f);
  ht = ht_alloc(meta.log2size, meta.size_key, meta.size_entry);
  fread(ht->entries, ht->size_entry, ht->size_proper, f);
  return ht;
}
