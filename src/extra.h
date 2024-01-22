
/**
   Global variables
 */

char dir_out[80]; // name of output dir
char file_param[80] = "params.in"; // default param filename


/** 
    Strings utilities
*/

#include <string.h>

char *strdup(const char *s) {
  size_t size = strlen(s) + 1;
  char *p = malloc(size);
  if (p) {
    memcpy(p, s, size);
  }
  return p;
}


void trim_whitespace(char *s) {
  const char *d = s;
  do {
    while (*d == ' ')
      ++d;
  } while ((*s++ = *d++));
  *s = '\0';
}


/** 
    Lists utilities
// from https://fr.wikibooks.org/wiki/Structures_de_donn%C3%A9es_en_C/Les_listes_simples
*/

typedef struct{
  char *name;
  void *data;
  char *type;
  int len; 
  
} List;


/* List * list_create (void *data, char *name, char *type) { */
/*   List *list = malloc(1 * sizeof *list); */
/*   list[0].data = data; */
/*   list[0].name = strdup(name); */
/*   list[0].type = strdup(type); */
/*   list[0].len = 1; */

/*   return list; */
/* } */


List * list_append(List *list, void *data, char *name, char *type) {


  if (list == NULL){ // creation of the list
    /* fprintf (stdout,"creation of the list\n"); */
    list = malloc(1 * sizeof *list);
    list[0].len = 0;
  }
  else{
    /* fprintf (stdout,"append to list\n");       */
    list = realloc(list, (list[0].len+1) * sizeof *list);
  }

  int len = list[0].len;
  list[len].data = data;
  list[len].name = strdup(name);
  list[len].type = strdup(type);
  list[0].len += 1; // store list length in list[0]

  return list;

}


void list_free(List *list){

  for (int i = 0; i < list[0].len; i++)
    free(list[i].name);
  
  free(list);
}


/**
   "Namelist" utilities
*/

void str2array(char *tmps2, double *array){

  // overwrite newline character with string terminator
  char *newline = strchr( tmps2, '\n' );
  if ( newline )
    *newline = 0;

  char *p;
  int n = 0;

  p = strtok(tmps2,"[,]");
  while (p != NULL){
    array[n] = atof(p);
    if (print){fprintf (stdout, "array[%d] = %g \n", n, array[n]);}
    p = strtok(NULL, ",");
    n += 1;
  }
}


void read_params(List *params, char *path2file)
{
  FILE *fp;
  if ((fp = fopen(path2file, "rt"))) {
    char tempbuff[300];
    while(fgets(tempbuff,300,fp)) { // loop over file lines
      trim_whitespace(tempbuff);
      char *tmps1 = strtok(tempbuff, "=");
      char *tmps2 = strtok(NULL, "=");

      //  ParamItem * d2 = params[0];
      for (int i = 0; i < params[0].len; i++) { // loop over parameters
        if (strcmp(params[i].name, tmps1) == 0) {
          //check for type and assign value
          if (strcmp(params[i].type, "int") == 0) {
            *( (int*) params[i].data) = atoi(tmps2);
            if (print){fprintf (stdout, "scan param %s: %d\n", params[i].name, *( (int*) params[i].data));}
            }
          else if (strcmp(params[i].type, "double") == 0){
            *( (double*) params[i].data) = atof(tmps2);
            if (print){fprintf (stdout, "scan param %s: %e\n", params[i].name, *( (double*) params[i].data));}
            }
          else if (strcmp(params[i].type, "array") == 0){
            if (print){fprintf (stdout, "scan param %s:\n", params[i].name);}
            str2array(tmps2, (double*) params[i].data);
          }
        }
      }

    }
    fclose(fp);
  } else {
    fprintf(stdout, "file %s not found\n", path2file);
    exit(0);
  }
}


/**
   Create output directory and copy input parameter file for backup
*/

// for mkdir
#include <sys/stat.h>
#include <sys/types.h>

void create_outdir()
{
  if (pid() == 0) {
    for (int i=1; i<10000; i++) {
      sprintf(dir_out, "outdir_%04d/", i);
      if (mkdir(dir_out, 0777) == 0) {
        fprintf(stdout,"Writing output in %s\n",dir_out);        
        break;
      }
    }
  }
#if _MPI
  MPI_Bcast(&dir_out, 80, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif
}


void backup_config(char *path2file)
{
  if (pid() == 0) {
    fprintf(stdout, "Backup config\n");
    char ch;
    char name[90];
    sprintf (name,"%sparams.in", dir_out);
    FILE *source = fopen(path2file, "r");
    FILE *target = fopen(name, "w");
    while ((ch = fgetc(source)) != EOF)
      fputc(ch, target);
    fclose(source);
    fclose(target);
  }
}


/**
   Copy file
   from: https://stackoverflow.com/questions/29079011/copy-file-function-in-c
*/
void backup_file(char *FileSource, char *dir_out)
{
  //  if (pid() == 0) {
  char FileDestination[100];
  sprintf (FileDestination,"%s%s", dir_out, FileSource);

  char    c[4096]; // or any other constant you like
  FILE    *stream_R = fopen(FileSource, "r");
  FILE    *stream_W = fopen(FileDestination, "w");   //create and write to file

  while (!feof(stream_R)) {
    size_t bytes = fread(c, 1, sizeof(c), stream_R);
    if (bytes) {
      fwrite(c, 1, bytes, stream_W);
    }
  }

  //close streams
  fclose(stream_R);
  fclose(stream_W);
  //  }
}
