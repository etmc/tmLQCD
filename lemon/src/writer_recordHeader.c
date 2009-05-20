#include <config.h>
#include <lemon.h>
#include "internal_writeRecordBinaryHeader.static"

int lemonWriteRecordHeader(LemonRecordHeader *props, LemonWriter* writer)
{
  int result;

  if ((writer == (LemonWriter *)NULL) || (props == (LemonRecordHeader *)NULL))
    return LEMON_ERR_PARAM;

  /* Make sure header is expected now */
  if (writer->header_nextP != 1)
    return LEMON_ERR_HEADER_NEXT;

  /* If last record ended a message, this one must begin a new one */
  /* If last record did not end a message, this one must not begin one */
  /* Since we allow appending to a file and we don't reread it to check
     the state of the last flag, we don't do this check for the first
     record written. */
  if (writer->first_record != 1 && writer->isLastP != props->MB_flag )
    return LEMON_ERR_MBME;

  result = writeLemonRecordBinaryHeader(writer, props);

  /* Set new writer state */
  writer->isLastP      = props->ME_flag;
  writer->first_record = props->MB_flag;
  writer->bytes_total  = props->data_length;
  writer->header_nextP = 0;

  return result;
}
