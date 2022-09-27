#ifndef _NEUTROOTHANDLERS_H_
#define _NEUTROOTHANDLERS_H_

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

extern "C"
{
  int rootopenfile_(char * /* filename */, char * /* option */, int , int);
  int rootwritefile_(char * /* filename */, int);
  int rootclosefile_(char * /* filename */, int);
};

class NeutRootHandlers
{
 public:

  int open(char *filename, char *opt);
  int write(char *filename);
  int close(char *filename);

  void nulltermstr(char *str, int len);

  TTree *maketree(char *filename, char *treename, char *title);
  TTree *attachtree(char *filename, char *treename);

};

#endif /* _NEUTROOTHANDLERS_H_ */
