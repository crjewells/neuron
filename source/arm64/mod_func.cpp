#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _double_synapse_reg(void);
extern void _gabaa_reg(void);
extern void _hh3k_reg(void);
extern void _hh3na_reg(void);
extern void _nmdaSyn_reg(void);
extern void _nmda_dsyn_reg(void);
extern void _syn_placeholder_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"double_synapse.mod\"");
    fprintf(stderr," \"gabaa.mod\"");
    fprintf(stderr," \"hh3k.mod\"");
    fprintf(stderr," \"hh3na.mod\"");
    fprintf(stderr," \"nmdaSyn.mod\"");
    fprintf(stderr," \"nmda_dsyn.mod\"");
    fprintf(stderr," \"syn_placeholder.mod\"");
    fprintf(stderr, "\n");
  }
  _double_synapse_reg();
  _gabaa_reg();
  _hh3k_reg();
  _hh3na_reg();
  _nmdaSyn_reg();
  _nmda_dsyn_reg();
  _syn_placeholder_reg();
}

#if defined(__cplusplus)
}
#endif
