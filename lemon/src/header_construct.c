#include <config.h>
#include <lemon.h>
#include <string.h>

LemonRecordHeader *lemonCreateHeader(int MB_flag, int ME_flag, char *type, uint64_t reclen)
{
  LemonRecordHeader *result;
  size_t type_length;

  result = (LemonRecordHeader*)malloc(sizeof(LemonRecordHeader));
  if (result == (LemonRecordHeader*)NULL)
    return NULL;

  type_length = strlen(type);
  result->type = (char*)malloc(sizeof(char) * (type_length + 1));
  if (result->type == (char*)NULL)
  {
    free(result);
    return NULL;
  }
  strcpy(result->type, type); /* Will be assumed to be null-terminated from now on */

  /* Fill out the rest of the fields */
  result->ME_flag = ME_flag;
  result->MB_flag = MB_flag;
  result->data_length = reclen;
  result->lemon_version = LEMON_VERSION;

  return result;
}
