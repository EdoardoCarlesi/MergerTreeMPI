#include <stdlib.h>

#ifndef PTREE_IO_H 
#define PTREE_IO_H
int      read_particles         (char filename[MAXSTRING], int isimu);

int      read_positions         (char filename[MAXSTRING], int isimu);

int      write_mtree            (int isimu0, char OutFile[MAXSTRING]);

int assign_input_files_to_tasks(char *partList, char *haloList, char *tempDir, char ***locPartFile, char ***locHaloFile, int nFiles);

int	 assign_output_files_names	(char *outList, char **outSuffix, int nFiles);
#endif
