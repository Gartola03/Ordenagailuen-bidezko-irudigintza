#include "obj.h"

void print_matrizea(char *str, double *m);

double modulua(double *v);

void biderkatu_matrizeak(double *emaitzaptr, double *eskerptr, double *eskuinptr);

double biderketa_eskalarra(double *v1, double *v2);

void matrizen_eta_bektorearen_biderketa(double *v, double *m_ptr, double *k);

void puntua_sortu(punto *p, vertex erpina);

void biderketa_bektoriala(double *v, double *v1, double *v2);

void bektorea_normalizatu(double *v);

void bektore_normala_kalkulatu(double *nptr, punto p1, punto p2, punto p3);

void bektore_normala_guztiak_kalkulatu(object3d *optr);

//matrize aldaketa

void erreferentzia_sistemaren_aldaketa(double *mesa, object3d *optr);

void modelview_kalkulatu(double *modelview, double *mesa, object3d *auxptr);

void kamara_arreta_puntuari_begiratu(object3d *selk_ptr, object3d *sel_ptr);

void biraketa_matrizea_analisi_moduan(object3d *selk_ptr, object3d *sel_ptr , double *mlag, char biraketa, double sinalfa, double cosalfa);

void perspektiba_proiekzioa(double *mperspektiba, double l, double r, double b, double t, double n, double f);

//marraztu

void gouraud_bektoreak_marraztu(object3d *optr, double *modelview, double *mperspektiba, int proiekzioa, int indexx);

void flat_bektoreak_marraztu(object3d *optr, double *modelview, double *mperspektiba, int proiekzioa, int i);


//puntuak

void mxp(punto *pptr, double m[16], punto p);

int mxpp(punto *pptr, double m[16], punto p);

void objektuari_aldaketa_sartu_ezk(double m[16], object3d *optr);

void objektuari_aldaketa_sartu_esk(double m[16], object3d *optr);

//argiak

void light_bektore_cam(punto p1, light *light_table, double *modelview);

void h_kalkulatu(double *H, double *vkamara, light *light_table);

void argi_kalkulua(color3 *I, object3d *optr, color3 Ia, int aurp, int gouraud, int ind, light *light_table, double *H);