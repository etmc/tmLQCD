/* Display contents of a LEMON formated file */
/* Balint Joo 2003 */
/* C. DeTar 10/26/04 reformatted output */
/* A. Deuzeman 21/04/09 MPI version */

#include <stdio.h>
#include <stdlib.h>
#include <lemon.h>
#include <mpi.h>
#define MAX_BYTES 0x1000

/* Scan for non-ASCII characters */
/* Return true if all characters are ASCII */
int all_ascii(char *buf, size_t length)
{
  size_t i;

  for(i = 0; i < length; i++)
    if(0x80 & buf[i])return 0;
    return 1;
}

int main(int argc, char *argv[])
{
  char* data_buf;

  LemonReader *reader;
  MPI_File fp;
  MPI_Info info;
  int status;
  uint64_t nbytes, read_bytes;
  int msg,rec,first;
  char *lemon_type;
  size_t bytes_pad;
  int MB_flag, ME_flag;

  MPI_Init(&argc, &argv);

  if( argc != 2 )
  {
    fprintf(stderr, "Usage: %s <lemon_file>\n", argv[0]);
    return EXIT_FAILURE;
  }

  MPI_Info_create(&info);
  MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, info, &fp);

  reader = lemonCreateReader(&fp, MPI_COMM_WORLD);
  if (reader == (LemonReader *)NULL)
  {
    fprintf(stderr, "Unable to open LemonReader.\n");
    return EXIT_FAILURE;
  }

  msg = 0; first = 1; rec = 0;
  while ((status = lemonReaderNextRecord(reader)) != LEMON_EOF)
  {
    if (status != LEMON_SUCCESS)
    {
      fprintf(stderr, "lemonReaderNextRecord returned status = %d.\n", status);
      return EXIT_FAILURE;
    }

    nbytes    = lemonReaderBytes(reader);
    lemon_type = lemonReaderType(reader);
    bytes_pad = lemonReaderPadBytes(reader);
    MB_flag   = lemonReaderMBFlag(reader);
    ME_flag   = lemonReaderMEFlag(reader);

    if (MB_flag == 1 || first)
    {
      first = 0;
      rec = 0;
      msg++;
    }

    rec++;

    printf("\n\n");
    printf("Message:        %d\n", msg);
    printf("Record:         %d\n", rec);
    printf("Type:           %s\n", lemon_type);
    printf("Data Length:    %llu\n", (unsigned long long)nbytes);
    printf("Padding Length: %lu\n", (unsigned long)bytes_pad);
    printf("MB flag:        %d\n", MB_flag);
    printf("ME flag:        %d\n", ME_flag);

    if (nbytes >= MAX_BYTES)
    {
      printf("Data:           [Long record skipped]\n");
      lemonReaderCloseRecord(reader);
      continue;
    }

    data_buf = (char *)malloc((size_t)nbytes+1);
    if (data_buf == (char *)NULL)
    {
      fprintf(stderr, "Couldn't malloc data buf.\n");
      return EXIT_FAILURE;
    }

    read_bytes = nbytes;
    status = lemonReaderReadData((void *)data_buf, &read_bytes, reader);

    if( status < 0 )
    {
      if( status != LEMON_EOR )
      {
        fprintf(stderr, "LEMON read error occurred: status= %d  %llu bytes wanted, %llu read.\n",
                status, (unsigned long long)nbytes, (unsigned long long)read_bytes);
        return EXIT_FAILURE;
      }
    }

    data_buf[nbytes]='\0';
    if(!all_ascii(data_buf, nbytes))
      printf("Data:           [Binary data].\n");
    else
      printf("Data:           \"%s\" \n", data_buf);

    free(data_buf);
    lemonReaderCloseRecord(reader);
  }

  lemonDestroyReader(reader);
  MPI_File_close(&fp);
  MPI_Finalize();

  return EXIT_SUCCESS;

}
