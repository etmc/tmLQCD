#include <config.h>

#include <stdio.h>

#include <lemon.h>

#include "internal_readAndParseHeader.static"
#include "internal_padding.static"

int lemonReaderNextRecord(LemonReader *reader)
{
  int err = 0;

  if (reader == (LemonReader *)NULL)
  {
    fprintf(stderr, "Node %d reports in lemonReaderNextRecord: got NULL pointer for reader argument.\n",
            reader->my_rank);
    return LEMON_ERR_PARAM;
  }

  if (reader->header_nextP != 1)
    lemonReaderCloseRecord(reader);

  err = readAndParseHeader(reader);

  if (err != LEMON_SUCCESS)
    return err;

  reader->is_last = reader->curr_header->ME_flag;
  reader->bytes_total = reader->curr_header->data_length;

  reader->bytes_pad = lemon_padding(reader->bytes_total);
  reader->header_nextP = 0;

  return LEMON_SUCCESS;
}
