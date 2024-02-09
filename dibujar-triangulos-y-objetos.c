//	Program developed by
//
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc dibujar-triangulos-y-objetos.c -lGL -lGLU -lglut -lm
//
//
//
#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "obj.h"
#include "funtzioak.h"

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy;

int indexx;
object3d *foptr;
object3d *sel_ptr;
object3d *fkptr;
object3d *selk_ptr;
light *light_table;
int num_lights;
color3 Ia;
color3 I;
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;
int kamara;
int objektua_sortu;
int ikuspegia;
int analisia;
int bektoreak;
int back_culling;
int proiekzioa;
int gouraud;
int argia_aukeratu;
double mesa[16];
double modelview[16];
double mperspektiba[16];

double H[12];

char fitxiz[100];


void info (){
    printf("\n");
    printf("INFORMAZIOA:\n");
    if (light_table[0].onoff == 1){
        printf("-EGUZKIA PIZTUTA dago.\n");
    } else{
        printf("-EGUZKIA ITZALITA dago.\n");
    }

    if (light_table[1].onoff == 1){
        printf("-BONBILLA PIZTUTA dago.\n");
    } else{
        printf("-BONBILLA ITZALITA dago.\n");
    }

    if (light_table[2].onoff == 1){
        printf("-OBJEKTUAREN BONBILLA PIZTUTA dago.\n");
    } else{
        printf("-OBJEKTUAREN BONBILLA ITZALITA dago.\n");
    }

    if (light_table[3].onoff == 1){
        printf("-KAMERAREN BONBILLA PIZTUTA dago.\n");
    } else{
        printf("-KAMERAREN BONBILLA ITZALITA dago.\n");
    }

    if (gouraud == 1){
        printf("-GOURAUD erabiltzen ari da.\n");
    } else{
        printf("-FLAT erabiltzen ari da.\n");
    }

    if (proiekzioa == 1){
        printf("-PROIEKZIOA PARALELOAN ikusten ari zara.\n");
    } else{
        printf("-PERSPEKTIBAN PROIEKZIOA ikusten ari zara.\n");
    }

    if (ikuspegia == 1){
        printf("-KAMERAREN IKUSPEGIAN ikusten ari zara.\n");
    } else{
        printf("-OBJEKTUAREN IKUSPEGIAN ikusten ari zara.\n");
    }

    if (kamara == 1){
        printf("-KAMERA AUKERATZEN ari zara.\n");
    } else if (kamara == 0){
        printf("-OBJEKTUA AUKERATZEN ari zara.\n");
    } else{
        printf("-ARGIA AUKERATZEN ari zara.\n");
        switch (argia_aukeratu){
        case 0:
            printf(":Eguzkia\n");
            break;
        case 1:
            printf(":Bonbilla\n");
            break;
        case 2:
            printf(":Objektuaren fokoa\n");
            break;
        case 3:
            printf(":Kameraren fokoa\n");
            break;
        }
    }

    if (ald_lokala == 1){
        printf("-ALDAKETA LOKALA (hegan moduan) eragiten ari zara.\n");
    } else{
        printf("-ALDAKETA GLOBALA (analisi moduan) eragiten ari zara.\n");
    }

    if (aldaketa == 't'){
        printf("-TRANSLAZIOAK eragingo dizkiozu.\n");
    } else{
        printf("-BIRAKETAK eragingo dizkiozu.\n");
    }

    if (bektoreak == 1){
        printf("-BEKTOREAK AKTIBATUTA daude.\n");
    } else{
        printf("-BEKTOREAK DESAKTIBATUTA daude.\n");
    }

    if (back_culling == 0){
        printf("-ATZE-AURPEGIEN MARRAZKETA egiten ari zara.\n");
    } else{
        printf("-ATZE-AURPEGIEN EZABAKETA egiten ari zara.\n");
    }
    printf("-KAMERAREN KOKAPENA: ");
    printf("%lf, %lf, %lf", selk_ptr->mptr->m[3], selk_ptr->mptr->m[7], selk_ptr->mptr->m[11]);
    printf("\n");
}

int ikusten_da(object3d *aukptr, object3d *optr, double *modelview, int i){
    punto p1;
    int j, k, ind1;
    double n[3];
    double vkamara[3];

    ind1 = optr->face_table[i].vertex_ind_table[0];
    puntua_sortu(&p1, optr->vertex_table[ind1]);
    mxp(&p1,modelview,p1);

    if (proiekzioa == 1){
        vkamara[0] = aukptr->mptr->m[2];
        vkamara[1] = aukptr->mptr->m[6];
        vkamara[2] = aukptr->mptr->m[10];

    }else{
        vkamara[0] = aukptr->mptr->m[3];
        vkamara[1] = aukptr->mptr->m[7];
        vkamara[2] = aukptr->mptr->m[11];

        vkamara[0] -= p1.x;
        vkamara[1] -= p1.y;
        vkamara[2] -= p1.z; 
    }
    
    light_bektore_cam(p1, light_table, &(mesa[0]));
    h_kalkulatu(&(H[0]), &(vkamara[0]), light_table);

    matrizen_eta_bektorearen_biderketa(&(optr->face_table[i].Ncam[0]),&(modelview[0]) ,&(optr->face_table[i].N[0]));
    bektorea_normalizatu(&(optr->face_table[i].Ncam[0]));
    if (optr != sel_ptr) glColor3ub(255, 255, 255);
    else glColor3ub(0, 255, 0);
    if (biderketa_eskalarra((&(optr->face_table[i].Ncam[0])), &(vkamara[0])) <= 0){
        if (back_culling == 1){
            return 1;
        }else{
            glColor3ub(255, 0, 0);
        }
    }
    return 0;
}

// TODO
// funtzio honek u eta v koordenatuei dagokien pointerra itzuli behar du.
// debe devolver el pointer correspondiente a las coordenadas u y v

unsigned char * color_textura(float u, float v)
{
    int desplazamendua,xind,yind;
    char * lag;

    xind = u*(dimx-1);
    xind = xind%dimx;
    yind = v*(dimy-1);
    yind = yind%dimy;
    yind = dimy-yind-1;
    desplazamendua = yind*dimx+xind;
    lag = (unsigned char *)bufferra;
    return(lag+3*desplazamendua);
}

// TODO
// lerroa marrazten du, baina testuraren kodea egokitu behar da
// dibuja una linea pero hay que codificar la textura

void  dibujar_linea_z(int linea,float c1x, float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v, color3 Iez, color3 Ies){
    float xkoord,zkoord;
    float u,v;
    float t,deltat;
    unsigned char r,g,b;
    unsigned char *colorv;

    glBegin( GL_POINTS );
    if(c1x != c2x){
        deltat = 1.0/(c2x-c1x);
        //Lerroan xkoord, zkoord, u eta v posizio berriak eguneratzen du. Zuzenaren koordenatu barizentrikoa erabiltzen du
        for (t = 1.0; t>=0;t-=deltat){
            xkoord = c1x*t +(1-t)*c2x;
            zkoord = c1z*t +(1-t)*c2z;
            u = c1u*t +(1-t)*c2u;
            v = c1v*t +(1-t)*c2v;
            //Texturaren posizioa jasotzen du
            /*
            colorv=  color_textura(u, v);
            r=colorv[0];
            g=colorv[1];
            b=colorv[2];
            */
            if (gouraud == 0){
                r = I.r;
                g = I.g;
                b = I.b;
            } else{
                r = Iez.r*t +(1-t)*Ies.r;
                g = Iez.g*t +(1-t)*Ies.g;
                b = Iez.b*t +(1-t)*Ies.b; 
                if (r > 255) r = 255;
                if (g > 255) g = 255;
                if (b > 255) b = 255;
            }
            // glColor3f [0-1]
            glColor3ub(r,g,b);
            glVertex3f(xkoord, linea, zkoord);
        }
    }
    glEnd();
}

//Triangeluan lerroan marrazteko behar duen bi ebakidura puntuak kalkulatzen du.
void ebaketa_kalkulatu(punto *Gptr, punto *Bptr, int h, punto *eptr, color3 *Ig, color3 *Ib, color3 *Iy)
{
    float alpha, beta;

    alpha = (float)(Gptr->y-h)/ (Gptr->y - Bptr->y);
    beta = 1.0-alpha;
    
    eptr -> y = h;
    eptr -> x = beta * Gptr->x + alpha * Bptr->x;
    eptr -> z = beta * Gptr->z + alpha * Bptr->z;
    eptr -> u = beta * Gptr->u + alpha * Bptr->u;
    eptr -> v = beta * Gptr->v + alpha * Bptr->v;
    
    Iy -> r = 1/(Gptr->y - Bptr->y) * (Ig->r * (eptr->y - Bptr->y) + Ib->r * (Gptr->y - eptr->y));
    Iy -> g = 1/(Gptr->y - Bptr->y) * (Ig->g * (eptr->y - Bptr->y) + Ib->g * (Gptr->y - eptr->y));
    Iy -> b = 1/(Gptr->y - Bptr->y) * (Ig->b * (eptr->y - Bptr->y) + Ib->b * (Gptr->y - eptr->y));

    if (Iy->r > 255) Iy->r = 255;
    if (Iy->g > 255) Iy->g = 255;
    if (Iy->b > 255) Iy->b = 255;
    
    if (Gptr->x > Bptr->x){
        eptr -> x = fmaxf(fminf(eptr->x, Gptr->x), Bptr->x);
    } else{
        eptr -> x = fmaxf(fminf(eptr->x, Bptr->x), Gptr->x);
    }
    
}

void dibujar_triangulo(object3d *optr, int aurpegia, int j){
    object3d *auxptr;
    punto *pgoiptr, *pbeheptr, *perdiptr, *pezkeptr, *peskuinptr, *plagptr;
    float x1,h1,z1,u1,v1,x2,h2,z2,u2,v2,x3,h3,z3,u3,v3;
    float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v;
    int linea, ind1, ind2, ind3;
    float cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
    punto p1,p2,p3;
    punto e1,e2,e;
    double n[3];
    int i,a,b,c;
    color3 I1,I2,I3;
    color3 *Igptr, *Ieptr, *Ibptr, *Ilagptr;
    color3 Ie1, Ie2;

    if (aurpegia >= optr->num_faces) return;
    
    ind1 = optr->face_table[aurpegia].vertex_ind_table[0];
    ind2 = optr->face_table[aurpegia].vertex_ind_table[j-1];
    ind3 = optr->face_table[aurpegia].vertex_ind_table[j];
    
    puntua_sortu(&p1, optr->vertex_table[ind1]);
    puntua_sortu(&p2, optr->vertex_table[ind2]);
    puntua_sortu(&p3, optr->vertex_table[ind3]);


    mxp(&p1,modelview,p1);
    mxp(&p2,modelview,p2);
    mxp(&p3,modelview,p3);

    
    if (gouraud == 1){
        argi_kalkulua(&I1, optr, Ia, aurpegia, gouraud, ind1, light_table, &(H[0]));
        argi_kalkulua(&I2, optr, Ia, aurpegia, gouraud, ind2, light_table, &(H[0]));
        argi_kalkulua(&I3, optr, Ia, aurpegia, gouraud, ind3, light_table, &(H[0]));
    } else{
        argi_kalkulua(&I,optr, Ia, aurpegia, gouraud, 0, light_table, &(H[0]));
    }

    if (proiekzioa == 0){
        a = mxpp(&p1,mperspektiba,p1);
        b = mxpp(&p2,mperspektiba,p2);
        c = mxpp(&p3,mperspektiba,p3);
        if (a == 1 || b == 1 || c == 1) return;
    }

    if (lineak == 1){
        glBegin(GL_POLYGON);
            glVertex3d(p1.x, p1.y, p1.z);
            glVertex3d(p2.x, p2.y, p2.z);
            glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        return;
    }

    if (p1.y > p2.y){
        pgoiptr = &p1;
        Igptr = &I1;
        pbeheptr = &p2;
        Ibptr = &I2;
    } else{
        pgoiptr = &p2;
        Igptr = &I2;
        pbeheptr = &p1;
        Ibptr = &I1;
    }
    if (p3.y > pgoiptr -> y){
        perdiptr = pgoiptr;
        Ieptr = Igptr;
        pgoiptr = &p3;
        Igptr = &I3;
    }else{
        if(p3.y < pbeheptr -> y){
            perdiptr = pbeheptr;
            Ieptr = Ibptr;
            pbeheptr = &p3;
            Ibptr = &I3;
        }
        else {
            perdiptr = &p3;
            Ieptr = &I3;
        }
    }

    if (pgoiptr->y != perdiptr->y){
        for (i = pgoiptr->y; i > perdiptr->y; i--) {
            ebaketa_kalkulatu(pgoiptr, perdiptr, i, &e1, Igptr, Ieptr, &Ie1);
            ebaketa_kalkulatu(pgoiptr, pbeheptr, i, &e2, Igptr, Ibptr, &Ie2);
            if (e1.x > e2.x)
                dibujar_linea_z(i, e2.x, e2.z, e2.u, e2.v, e1.x, e1.z, e1.u, e1.v, Ie2, Ie1);
            else
                dibujar_linea_z(i, e1.x, e1.z, e1.u, e1.v, e2.x, e2.z, e2.u, e2.v, Ie1, Ie2);
        } 
    }

    if (perdiptr->y != pbeheptr->y){
        for (i = perdiptr->y - 1; i > pbeheptr->y; i--) {
            ebaketa_kalkulatu(perdiptr, pbeheptr, i, &e1, Ieptr, Ibptr, &Ie1);
            ebaketa_kalkulatu(pgoiptr, pbeheptr, i, &e2, Igptr, Ibptr, &Ie2);
            if (e1.x > e2.x)
                dibujar_linea_z(i, e2.x, e2.z, e2.u, e2.v, e1.x, e1.z, e1.u, e1.v, Ie2, Ie1);
            else
                dibujar_linea_z(i, e1.x, e1.z, e1.u, e1.v, e2.x, e2.z, e2.u, e2.v, Ie1, Ie2);
        }
    }
}

static void marraztu(void){
    float u,v;
    int i,j;
    object3d *auxptr, *aukptr;

    // marrazteko objektuak behar dira
    // no se puede dibujar sin objetos
    //if (foptr == 0) return;

    // clear viewport...
    if (objektuak == 1) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
        else
        {
        if (denak == 0) glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
        }

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    //fokoa
    light_table[2].pos[0] = sel_ptr->mptr->m[3];
    light_table[2].pos[1] = sel_ptr->mptr->m[7];
    light_table[2].pos[2] = sel_ptr->mptr->m[11];
    light_table[2].foko[0] = -sel_ptr->mptr->m[2];
    light_table[2].foko[1] = -sel_ptr->mptr->m[6];
    light_table[2].foko[2] = -sel_ptr->mptr->m[10];

    light_table[3].foko[0] = -selk_ptr->mptr->m[2];
    light_table[3].foko[1] = -selk_ptr->mptr->m[6];
    light_table[3].foko[2] = -selk_ptr->mptr->m[10];
    bektorea_normalizatu(&(light_table[3].foko[0]));
    

    if (proiekzioa == 1){
        glOrtho(-500.0, 500.0, -500.0, 500.0, 0.0, 500.0);
    } else{
        glOrtho(-500.0, 500.0, -500.0, 500.0, -500.0, 500.0);
    }

    if (ikuspegia == 1){
        aukptr = selk_ptr;
    }
    else{
        aukptr = sel_ptr;
    }

    erreferentzia_sistemaren_aldaketa(&(mesa[0]),aukptr);

    //kamera
    if (fkptr != 0){
        if (objektuak == 1){
            if (denak == 1){
                for (auxptr =fkptr; auxptr != 0; auxptr = auxptr->hptr){
                    modelview_kalkulatu(&(modelview[0]), &(mesa[0]), auxptr);
                    if (gouraud == 1){
                        for (i = 0;  i < auxptr->num_vertices; i++){
                            matrizen_eta_bektorearen_biderketa(&(auxptr->vertex_table[i].Ncam[0]), &(modelview[0]), &(auxptr->vertex_table[i].N[0]));
                            bektorea_normalizatu(&(auxptr->vertex_table[i].Ncam[0]));
                        }
                    }
                    matrizen_eta_bektorearen_biderketa(&(light_table[2].foko[0]),&(mesa[0]),&(light_table[2].foko[0]));
                    bektorea_normalizatu(&(light_table[2].foko[0]));
                    for (i = 0; i < auxptr->num_faces; i++){
                        if (ikusten_da(aukptr, auxptr, &(modelview[0]),i) == 0) {
                            for (j=2; j<auxptr->face_table[i].num_vertices; j++){
                                dibujar_triangulo(auxptr,i, j);
                            }
                            if (lineak == 1 && bektoreak == 1 && gouraud == 0){
                                flat_bektoreak_marraztu(auxptr, &(mesa[0]), &(mperspektiba[0]), proiekzioa, i); 
                            }                  
                        }
                    }
                    if (lineak == 1 && bektoreak == 1 && gouraud == 1){
                        gouraud_bektoreak_marraztu(auxptr, &(mesa[0]), &(mperspektiba[0]), proiekzioa, -1);
                    }
                }
            } else{
                modelview_kalkulatu(&(modelview[0]), &(mesa[0]), selk_ptr);
                if (gouraud == 1){
                    for (i = 0;  i < selk_ptr->num_vertices; i++){
                        matrizen_eta_bektorearen_biderketa(&(selk_ptr->vertex_table[i].Ncam[0]), &(modelview[0]), &(selk_ptr->vertex_table[i].N[0]));
                        bektorea_normalizatu(&(selk_ptr->vertex_table[i].Ncam[0]));
                    }
                }
                matrizen_eta_bektorearen_biderketa(&(light_table[2].foko[0]),&(mesa[0]),&(light_table[2].foko[0]));
                bektorea_normalizatu(&(light_table[2].foko[0]));
                for (i =0; i < selk_ptr->num_faces; i++){
                    if (ikusten_da(aukptr, selk_ptr, &(modelview[0]), i) == 0) {
                        for (j=2; j<selk_ptr->face_table[i].num_vertices; j++){
                            dibujar_triangulo(selk_ptr,i, j);
                        }
                        if (lineak == 1 && bektoreak == 1 && gouraud == 0){
                            flat_bektoreak_marraztu(selk_ptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa,i); 
                        } 
                    }
                }
                if (lineak == 1 && bektoreak == 1 && gouraud == 1){
                    gouraud_bektoreak_marraztu(selk_ptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, -1);
                }
            }
        } else{
            modelview_kalkulatu(&(modelview[0]), &(mesa[0]), selk_ptr);
            if (gouraud == 1){
                for (i = 0;  i < selk_ptr->num_vertices; i++){
                    matrizen_eta_bektorearen_biderketa(&(selk_ptr->vertex_table[i].Ncam[0]), &(modelview[0]), &(selk_ptr->vertex_table[i].N[0]));
                    bektorea_normalizatu(&(selk_ptr->vertex_table[i].Ncam[0]));
                }
            }
            matrizen_eta_bektorearen_biderketa(&(light_table[2].foko[0]),&(mesa[0]),&(light_table[2].foko[0]));
            bektorea_normalizatu(&(light_table[2].foko[0]));
            if (ikusten_da(aukptr, selk_ptr, &(modelview[0]), indexx) == 0) {
                dibujar_triangulo(selk_ptr,indexx, 2);
            }
            if (lineak == 1 && bektoreak == 1 && gouraud == 0){
                flat_bektoreak_marraztu(selk_ptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, indexx); 
            }
            else{
                gouraud_bektoreak_marraztu(selk_ptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, indexx);
                    
            }
        }
    }

    //objektua
    if (foptr != 0){
        if (objektuak == 1){
            if (denak == 1){
                for (auxptr = foptr; auxptr != 0; auxptr = auxptr->hptr){
                    modelview_kalkulatu(&(modelview[0]), &(mesa[0]), auxptr);
                    if (gouraud == 1){
                        for (i = 0;  i < auxptr->num_vertices; i++){
                            matrizen_eta_bektorearen_biderketa(&(auxptr->vertex_table[i].Ncam[0]), &(modelview[0]), &(auxptr->vertex_table[i].N[0]));
                            bektorea_normalizatu(&(auxptr->vertex_table[i].Ncam[0]));
                        }
                    }
                    matrizen_eta_bektorearen_biderketa(&(light_table[2].foko[0]),&(mesa[0]),&(light_table[2].foko[0]));
                    bektorea_normalizatu(&(light_table[2].foko[0]));
                    for (i =0; i < auxptr->num_faces; i++){
                        if (ikusten_da(aukptr, auxptr, &(modelview[0]), i) == 0) {
                            for (j=2; j<auxptr->face_table[i].num_vertices; j++){
                                dibujar_triangulo(auxptr,i, j);
                            }
                            if (lineak == 1 && bektoreak == 1 && gouraud == 0){
                                flat_bektoreak_marraztu(auxptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, i); 
                            }
                        }
                    }
                    if (lineak == 1 && bektoreak == 1 && gouraud == 1){
                        gouraud_bektoreak_marraztu(auxptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, -1);
                    }
                }
            } else{
                modelview_kalkulatu(&(modelview[0]), &(mesa[0]), sel_ptr);
                if (gouraud == 1){
                    for (i = 0;  i < sel_ptr->num_vertices; i++){
                        matrizen_eta_bektorearen_biderketa(&(sel_ptr->vertex_table[i].Ncam[0]), &(modelview[0]), &(sel_ptr->vertex_table[i].N[0]));
                        bektorea_normalizatu(&(sel_ptr->vertex_table[i].Ncam[0]));
                    }
                }
                matrizen_eta_bektorearen_biderketa(&(light_table[2].foko[0]),&(mesa[0]),&(light_table[2].foko[0]));
                bektorea_normalizatu(&(light_table[2].foko[0]));
                for (i =0; i < sel_ptr->num_faces; i++){
                    if (ikusten_da(aukptr, sel_ptr, &(modelview[0]), i) == 0) {
                        for (j=2; j<auxptr->face_table[i].num_vertices; j++){
                            dibujar_triangulo(auxptr,i, j);
                        }
                        if (lineak == 1 && bektoreak == 1 && gouraud == 0){
                            flat_bektoreak_marraztu(auxptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, i); 
                        }
                    }
                }
                if (lineak == 1 && bektoreak == 1 && gouraud == 1){
                    gouraud_bektoreak_marraztu(sel_ptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, -1);
                }
            }
        } else{
            modelview_kalkulatu(&(modelview[0]), &(mesa[0]), sel_ptr);
            if (gouraud == 1){
                for (i = 0;  i < auxptr->num_vertices; i++){
                    matrizen_eta_bektorearen_biderketa(&(auxptr->vertex_table[i].Ncam[0]), &(modelview[0]), &(auxptr->vertex_table[i].N[0]));
                    bektorea_normalizatu(&(auxptr->vertex_table[i].Ncam[0]));
                }
            }
            matrizen_eta_bektorearen_biderketa(&(light_table[2].foko[0]),&(mesa[0]),&(light_table[2].foko[0]));
            bektorea_normalizatu(&(light_table[2].foko[0]));
            if (ikusten_da(aukptr, sel_ptr, &(modelview[0]), indexx) == 0){
                dibujar_triangulo(sel_ptr,indexx, 2);
            }
            if (lineak == 1 && bektoreak == 1 && gouraud == 0){
                flat_bektoreak_marraztu(auxptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, indexx); 
            } else{
                gouraud_bektoreak_marraztu(sel_ptr, &(modelview[0]), &(mperspektiba[0]), proiekzioa, indexx);
            }
        }
    }
    glFlush();
    info(); 
}

void read_from_file(char *fitx){
    int i,n,retval;
    object3d *optr;

    printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (object3d *)malloc(sizeof(object3d));
    retval = read_wavefront(fitx, optr);

    if (retval !=0) {
        printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n",fitxiz);
        free(optr);
    } else {
        printf("objektuaren matrizea...\n");
        optr->mptr = (mlist *)malloc(sizeof(mlist));
        for (i=0; i<16; i++) optr->mptr->m[i] =0.0;
        optr->mptr->m[0] = 1.0;
        optr->mptr->m[5] = 1.0;
        optr->mptr->m[10] = 1.0;
        optr->mptr->m[15] = 1.0;
        optr->mptr->hptr = 0;

        bektore_normala_guztiak_kalkulatu(optr);
        
        printf("objektu zerrendara doa informazioa...\n");
        if (objektua_sortu == 1){
            optr->hptr = fkptr;
            fkptr = optr;
            selk_ptr = optr;
        }else{
            optr->hptr = foptr;
            foptr = optr;
            sel_ptr = optr;
        }
    }
    printf("datuak irakurrita\nLecura finalizada\n");
}

void argiak_sortu(){
    light_table = (light *) malloc(num_lights * sizeof (light));

    //onoff
    light_table[0].onoff = 0;
    light_table[1].onoff = 0;
    light_table[2].onoff = 0;
    light_table[3].onoff = 0;

    //type
    light_table[0].type = 0;
    light_table[1].type = 1;
    light_table[2].type = 2;
    light_table[3].type = 3;

    //nahi dudan koloreak
    light_table[0].I.r = 200;
    light_table[0].I.g = 200;
    light_table[0].I.b = 200;

    light_table[1].I.r = 250;
    light_table[1].I.g = 50;
    light_table[1].I.b = 150;

    light_table[2].I.r = 75;
    light_table[2].I.g = 75;
    light_table[2].I.b = 75;

    light_table[3].I.r = 155; 
    light_table[3].I.g = 120;
    light_table[3].I.b = 75;

    //nahi dudan posizioa
    light_table[0].pos[0] = 400;
    light_table[0].pos[1] = 400;
    light_table[0].pos[2] = 400;

    light_table[1].pos[0] = 0;
    light_table[1].pos[1] = 200;
    light_table[1].pos[2] = 0;

    light_table[2].pos[0] = sel_ptr->mptr->m[3];
    light_table[2].pos[1] = sel_ptr->mptr->m[7];
    light_table[2].pos[2] = sel_ptr->mptr->m[11];

    light_table[3].pos[0] = selk_ptr->mptr->m[3];
    light_table[3].pos[1] = selk_ptr->mptr->m[7];
    light_table[3].pos[2] = selk_ptr->mptr->m[11];

    //direkzioa
    light_table[0].dir[0] = 0.1;
    light_table[0].dir[1] = 0.2;
    light_table[0].dir[2] = 0.7;

    light_table[1].dir[0] = 0;
    light_table[1].dir[1] = -1;
    light_table[1].dir[2] = 0;

    //irerkiera
    light_table[0].aperture = 0;
    light_table[1].aperture = 0;
    light_table[2].aperture = 3.14/7;
    light_table[3].aperture = 3.14/6;//30

}

void x_aldaketa(int dir){

    double p, alfa, sinalfa, cosalfa;
    object3d *auxptr;
    char biraketa;
    int i;
    double mlag[16];
    for (i=0; i<16; i++) mlag[i] = 0.0;

    mlag[0] = 1.0;
    mlag[5] = 1.0;
    mlag[10] = 1.0;
    mlag[15] = 1.0;

    p=2.0;
    alfa=0.1;
    sinalfa = sin(alfa);
    cosalfa = cos(alfa);
    if (dir == 0){
        p=-p;
        alfa = -alfa;
        sinalfa = sin(alfa);
        cosalfa = cos(alfa);
    }
    if (kamara == 2){
        switch (argia_aukeratu){
        case 0:
            mlag[5]=cosalfa;
            mlag[6]=-sinalfa;
            mlag[9]=sinalfa;
            mlag[10]=cosalfa;
            matrizen_eta_bektorearen_biderketa(&(light_table[0].pos[0]), &(mlag[0]), &(light_table[0].pos[0]));
            matrizen_eta_bektorearen_biderketa(&(light_table[0].dir[0]), &(mlag[0]), &(light_table[0].dir[0]));
            break;
        case 1:
            light_table[1].pos[0]+=p;
            break;
        case 2:
            light_table[2].aperture+=alfa;
            if (light_table[2].aperture > 1.2211) light_table[2].aperture = 1.2211;
            if (light_table[2].aperture < 0) light_table[2].aperture = 0;
            break;
        case 3:
            light_table[3].aperture+=alfa;
            if (light_table[3].aperture > 1.2211) light_table[3].aperture = 1.2211;
            if (light_table[3].aperture < 0) light_table[3].aperture = 0;
            break;
        } 
    } else{
        if (kamara == 1){
            auxptr = selk_ptr;
            mlag[0]=cosalfa;
            mlag[2]=sinalfa;
            mlag[8]=-sinalfa;
            mlag[10]=cosalfa;
            if (analisia == 1){
                biraketa = 'x';
                biraketa_matrizea_analisi_moduan(selk_ptr, sel_ptr, &(mlag[0]), biraketa, sinalfa, cosalfa);
            }
        }else{
            auxptr = sel_ptr;
            mlag[3] = p;
            if(aldaketa=='r'){
                mlag[3]=0.0;
                mlag[5]=cosalfa;
                mlag[6]=-sinalfa;
                mlag[9]=sinalfa;
                mlag[10]=cosalfa;
            }
        }
        if(ald_lokala==1){
            objektuari_aldaketa_sartu_esk(mlag, auxptr);
        }else{
            objektuari_aldaketa_sartu_ezk(mlag, auxptr);
        }
    }
}


void y_aldaketa(int dir){
    double p, alfa, sinalfa, cosalfa;
    object3d *auxptr;
    char biraketa;
    int i;
    double mlag[16];
    for (i=0; i<16; i++) mlag[i] = 0.0;

    mlag[0] = 1.0;
    mlag[5] = 1.0;
    mlag[10] = 1.0;
    mlag[15] = 1.0;

    p=2.0;
    alfa=0.1;
    sinalfa = sin(alfa);
    cosalfa = cos(alfa);
    if (dir == 0){
        p=-p;
        alfa = -alfa;
        sinalfa = sin(alfa);
        cosalfa = cos(alfa);
    }
    if (kamara == 2){
        switch (argia_aukeratu){
        case 0:
            mlag[0]=cosalfa;
            mlag[2]=sinalfa;
            mlag[8]=-sinalfa;
            mlag[10]=cosalfa;
            matrizen_eta_bektorearen_biderketa(&(light_table[0].pos[0]), &(mlag[0]), &(light_table[0].pos[0]));
            matrizen_eta_bektorearen_biderketa(&(light_table[0].dir[0]), &(mlag[0]), &(light_table[0].dir[0]));
            break;
        case 1:
            light_table[1].pos[1]+=p;
            break;
        case 2:
            light_table[2].aperture+=alfa;
            if (light_table[2].aperture > 1.2211) light_table[2].aperture = 1.2211;
            if (light_table[2].aperture < 0) light_table[2].aperture = 0;
            break;
        case 3:
            light_table[3].aperture+=alfa;
            if (light_table[3].aperture > 1.2211) light_table[3].aperture = 1.2211;
            if (light_table[3].aperture < 0) light_table[3].aperture = 0;
            break;
        }
    } else{
        if (kamara == 1){
            auxptr = selk_ptr;
            mlag[5]=cosalfa;
            mlag[6]=-sinalfa;
            mlag[9]=sinalfa;
            mlag[10]=cosalfa;
            if (analisia == 1){
                biraketa = 'y';
                biraketa_matrizea_analisi_moduan(selk_ptr, sel_ptr, &(mlag[0]), biraketa, sinalfa, cosalfa);
            }
        }else{
            auxptr = sel_ptr;
            mlag[7] = p;
            if(aldaketa=='r'){
                mlag[0]=cosalfa;
                mlag[2]=sinalfa;
                mlag[7]=0.0;
                mlag[8]=-sinalfa;
                mlag[10]=cosalfa;
            }
        }
        if(ald_lokala==1){
            objektuari_aldaketa_sartu_esk(mlag, auxptr);
        }else{
            objektuari_aldaketa_sartu_ezk(mlag, auxptr);
        }
    }
}

void z_aldaketa(int dir){
    double p, alfa, sinalfa, cosalfa;
    object3d *auxptr;
    char biraketa;
    int i;
    double mlag[16];
    for (i=0; i<16; i++) mlag[i] = 0.0;

    mlag[0] = 1.0;
    mlag[5] = 1.0;
    mlag[10] = 1.0;
    mlag[15] = 1.0;

    p=2.0;
    alfa=0.1;
    sinalfa = sin(alfa);
    cosalfa = cos(alfa);
    auxptr = sel_ptr;
    if (kamara == 1){
        auxptr = selk_ptr;
    }

    if (dir == 0){
        p=-p;
        alfa = -alfa;
        sinalfa = sin(alfa);
        cosalfa = cos(alfa);
    }
    if (kamara == 2){
        switch (argia_aukeratu){
        case 0:
            mlag[0]=cosalfa;
            mlag[1]=-sinalfa;
            mlag[4]=sinalfa;
            mlag[5]=cosalfa;
            matrizen_eta_bektorearen_biderketa(&(light_table[0].pos[0]), &(mlag[0]), &(light_table[0].pos[0]));
            matrizen_eta_bektorearen_biderketa(&(light_table[0].dir[0]), &(mlag[0]), &(light_table[0].dir[0]));
            break;
        case 1:
            light_table[1].pos[2]+=p;
            break;
        case 2:
            light_table[2].aperture+=alfa;
            if (light_table[2].aperture > 1.2211) light_table[2].aperture = 1.2211;
            if (light_table[2].aperture < 0) light_table[2].aperture = 0;
            break;
        case 3:
            light_table[3].aperture+=alfa;
            if (light_table[3].aperture > 1.2211) light_table[3].aperture = 1.2211;
            if (light_table[3].aperture < 0) light_table[3].aperture = 0;
            break;
        }
    } else{
        mlag[11] = p;
        if(aldaketa=='r'){
            mlag[0]=cosalfa;
            mlag[1]=-sinalfa;
            mlag[4]=sinalfa;
            mlag[5]=cosalfa;
            mlag[11]=0.0;
            if (kamara == 1 && analisia == 1){
                biraketa = 'z';
                biraketa_matrizea_analisi_moduan(selk_ptr, sel_ptr, &(mlag[0]), biraketa, sinalfa, cosalfa);
            }
        }

        if (ald_lokala==1){
            objektuari_aldaketa_sartu_esk(mlag, auxptr);
        } else{
            if (aldaketa != 'r'){
                objektuari_aldaketa_sartu_esk(mlag, auxptr);
            } else{
                objektuari_aldaketa_sartu_ezk(mlag, auxptr);
            }
        }
    }
}

void eskala_aldaketa(int dir){
    double p;
    object3d *auxptr;
    int i;
    double mlag[16];
    for (i=0; i<16; i++) mlag[i] = 0.0;

    p = 1.2;
    if (dir == 0){
        p=1/p;
    }
    auxptr = sel_ptr;
    if (kamara == 1){
        auxptr = selk_ptr;
    }
    mlag[0] = p;
    mlag[5] = p;
    mlag[10] = p;
    mlag[15] = 1.0;


    if(ald_lokala==1){
        objektuari_aldaketa_sartu_esk(mlag, auxptr);
    } else{
        objektuari_aldaketa_sartu_ezk(mlag, auxptr);
    }

}

void undo(){
    object3d *auxptr;
    auxptr = sel_ptr;
    if (kamara == 1){
        auxptr = selk_ptr;
    }
    if(auxptr->mptr->hptr!=0){
        auxptr->mptr= auxptr->mptr->hptr;
    }
}

// This function will be called whenever the user pushes one key
static void teklatua (unsigned char key, int x, int y){
    int retval;
    int i;
    double alfa;
    FILE *obj_file;

    switch(key){
        /* aldatu???????????????????
        case 13:
            if (foptr != 0)  // objekturik ez badago ezer ez du egin behar
                            // si no hay objeto que no haga nada
                {
                indexx ++;  // azkena bada lehenengoa bihurtu
                        // pero si es el último? hay que controlarlo!
            if (indexx == sel_ptr->num_triangles)
                {
                indexx = 0;
                if ((denak == 1) && (objektuak == 0))
                    {
                    glClear( GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT );
                    glFlush();
                    }
                }
            }
            break;
        */
        case 'd':
            if (denak == 1) denak = 0;
                else denak = 1;
            break;
        case 'o':
            if (objektuak == 1) objektuak = 0;
                else objektuak = 1;
            break;
        case 'l':
            if (lineak == 1) lineak = 0;
                else lineak = 1;
            break;
        case 't':
                aldaketa = 't';
            break;
        case 'r':
            aldaketa = 'r';
            break;
        case 'g':
            if (ald_lokala == 1) {
                ald_lokala = 0; //aldaketa globala
                analisia = 1;
                if (kamara == 1){
                    kamara_arreta_puntuari_begiratu(selk_ptr, sel_ptr);
                }
            }
            else {
                ald_lokala = 1;
                analisia = 0;
            }
            break;
        case 'c':
            if (kamara == 0) {
                kamara = 1;
                if (analisia == 1){
                    kamara_arreta_puntuari_begiratu(selk_ptr, sel_ptr);
                }
            }else if (kamara == 1){
                //argiak
                kamara = 2;
            }else{
                kamara = 0;
            }
            break;
        case 'C':
            if (ikuspegia == 1) {
                ikuspegia = 0;
            }else {
                ikuspegia= 1;
            }
            break;
        case 'n':
            if (bektoreak == 1){
                bektoreak = 0;
            }else {
                bektoreak = 1;
            }
            break;
        case 'b':
            if (back_culling == 1){
                back_culling = 0;
            }
            else {
                back_culling= 1;
            }
            break;
        case 'p':
            if (proiekzioa == 1){
                proiekzioa = 0;
            }
            else {
                proiekzioa = 1;
            }
            break;
        case 'x':
            x_aldaketa(1);
            break;
        case 'y':
            y_aldaketa(1);
            break;
        case 'z':
            z_aldaketa(1);
            break;
        case 'X':
            x_aldaketa(0);
            break;
        case 'Y':
            y_aldaketa(0);
            break;
        case 'Z':
            z_aldaketa(0);
            break;
        case 's':
            eskala_aldaketa(0);
            break;
        case 'S':
            eskala_aldaketa(1);
            break;
        case 'u':
            undo();
            break;

        case 'F':
            if (gouraud == 1){
                gouraud = 0;
            }else {
                gouraud = 1;
            }
            break;
        case '0':
            if (light_table[0].onoff  == 1){
                light_table[0].onoff  = 0;
            }else {
                light_table[0].onoff  = 1;
            }
            break;

        case '1':
            if (light_table[1].onoff  == 1){
                light_table[1].onoff = 0;
            }else {
                light_table[1].onoff = 1;
            }
            break;

        case '2':
            if (light_table[2].onoff == 1){
                light_table[2].onoff = 0;
            }else {
                light_table[2].onoff = 1;
            }
            break;

        case '3':
            if (light_table[3].onoff == 1){
                light_table[3].onoff = 0;
            }else {
                light_table[3].onoff = 1;
            }
            break;

        case '+':
            alfa =0.1;
            if (kamara == 2){
                switch (argia_aukeratu){
                case 2:
                    light_table[2].aperture+=alfa;
                    if (light_table[2].aperture > 1.2211) light_table[2].aperture = 1.2211;
                    if (light_table[2].aperture < 0) light_table[2].aperture = 0;
                    break;
                case 3:
                    light_table[3].aperture+=alfa;
                    if (light_table[3].aperture > 1.2211) light_table[3].aperture = 1.2211;
                    if (light_table[3].aperture < 0) light_table[3].aperture = 0;
                    break;
                default:
                    break;
                }
            }
            break;

        case '-':
            alfa =-0.1;
            if (kamara == 2){
                switch (argia_aukeratu){
                case 2:
                    light_table[2].aperture+=alfa;
                    if (light_table[2].aperture > 1.2211) light_table[2].aperture = 1.2211;
                    if (light_table[2].aperture < 0) light_table[2].aperture = 0;
                    break;
                case 3:
                    light_table[3].aperture+=alfa;
                    if (light_table[3].aperture > 1.2211) light_table[3].aperture = 1.2211;
                    if (light_table[3].aperture < 0) light_table[3].aperture = 0;
                    break;
                default:
                    break;
                }
            }
            break;
        case 'f':
                /*Ask for file*/
                printf("\n");
                if (kamara != 2){
                    if(kamara == 1){
                        objektua_sortu = 1;
                        printf("Kamara bat sartuko duzu.\n");
                        printf("idatzi fitxategi izena:\n");
                        scanf("%s", &(fitxiz[0]));
                        read_from_file(fitxiz);
                        selk_ptr->rgb.b = 140.0;
                        selk_ptr->rgb.b = 45.0;
                        selk_ptr->rgb.b = 255.0;

                        selk_ptr->Ka.r = 0.9;
                        selk_ptr->Ka.g = 0.3;
                        selk_ptr->Ka.b = 0.2;
                        
                        selk_ptr->kd.r = 0.5922;
                        selk_ptr->kd.g = 0.0166;
                        selk_ptr->kd.b = 0.0;

                        selk_ptr->ks.r = 0.5974;
                        selk_ptr->ks.g = 0.2084;
                        selk_ptr->ks.b = 0.2084;

                        selk_ptr->ns = 100.2237;
                        indexx = 0;
                    }else{
                        objektua_sortu = 0;
                        printf("Objektu bat sartuko duzu.\n");
                        printf("idatzi fitxategi izena:\n");
                        scanf("%s", &(fitxiz[0]));
                        read_from_file(fitxiz);
                        sel_ptr->rgb.b = 140.0;
                        sel_ptr->rgb.b = 45.0;
                        sel_ptr->rgb.b = 255.0;

                        sel_ptr->Ka.r = 0.9;
                        sel_ptr->Ka.g = 0.3;
                        sel_ptr->Ka.b = 0.2;
                        
                        sel_ptr->kd.r = 0.5922;
                        sel_ptr->kd.g = 0.0166;
                        sel_ptr->kd.b = 0.0;

                        sel_ptr->ks.r = 0.5974;
                        sel_ptr->ks.g = 0.2084;
                        sel_ptr->ks.b = 0.2084;

                        sel_ptr->ns = 100.2237;
                        indexx = 0;
                    }

                } else{
                    printf("Ezin da argi bat sartu!");
                }
                break;
        /* case 'S':  // save to file
                printf("idatzi fitxategi izena\n");
                scanf("%s", &(fitxiz[0]));
                    if ((obj_file = fopen(fitxiz, "w")) == NULL)
                            {
                            printf("ezin fitxategia ireki\n");
                            }
                        else
                            {
                            for (i =0; i < sel_ptr->num_triangles; i++)
                                {
                                fprintf(obj_file,"t %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                                    sel_ptr->triptr[i].p1.x-250, sel_ptr->triptr[i].p1.y-250, sel_ptr->triptr[i].p1.z,
                                    sel_ptr->triptr[i].p1.u, sel_ptr->triptr[i].p1.v,
                                    sel_ptr->triptr[i].p2.x-250, sel_ptr->triptr[i].p2.y-250, sel_ptr->triptr[i].p2.z,
                                    sel_ptr->triptr[i].p2.u, sel_ptr->triptr[i].p2.v,
                                    sel_ptr->triptr[i].p3.x-250, sel_ptr->triptr[i].p3.y-250, sel_ptr->triptr[i].p3.z,
                                    sel_ptr->triptr[i].p3.u, sel_ptr->triptr[i].p3.v );
                                }
                            fclose(obj_file);
                            }
                    break; */
        case 9: /* <TAB> */

            if (foptr != 0) // objekturik gabe ez du ezer egin behar
                            // si no hay objeto no hace nada
                {
                if (kamara == 2){
                    argia_aukeratu++;
                    if (argia_aukeratu == 4) argia_aukeratu = 0;
                } else{
                    if (kamara == 1){
                        selk_ptr = selk_ptr->hptr;
                        /*The selection is circular, thus if we move out of the list we go back to the first element*/
                        if (selk_ptr == 0) selk_ptr = fkptr;

                        light_table[3].pos[0] = selk_ptr->mptr->m[3];
                        light_table[3].pos[1] = selk_ptr->mptr->m[7];
                        light_table[3].pos[2] = selk_ptr->mptr->m[11];

                        indexx =0; // the selected polygon is the first one

                    }else{
                        sel_ptr = sel_ptr->hptr;
                        /*The selection is circular, thus if we move out of the list we go back to the first element*/
                        if (sel_ptr == 0) sel_ptr = foptr;

                        light_table[2].pos[0] = sel_ptr->mptr->m[3];
                        light_table[2].pos[1] = sel_ptr->mptr->m[7];
                        light_table[2].pos[2] = sel_ptr->mptr->m[11];
                        indexx =0; // the selected polygon is the first one
                    }
                }
            }
            break;
        case 27:  // <ESC>
            exit( 0 );
            break;
        default:
            printf("%d %c\n", key, key );
        }

    // The screen must be drawn to show the new triangle
    glutPostRedisplay();
}

int main(int argc, char** argv){
    int retval;

	printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
	printf("Press <ESC> to finish\n");
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_RGB|GLUT_DEPTH );
	glutInitWindowSize ( 500, 500 );
	glutInitWindowPosition ( 100, 100 );
	glutCreateWindow( "KBG/GO praktika" );

	glutDisplayFunc( marraztu );
	glutKeyboardFunc( teklatua );
	/* we put the information of the texture in the buffer pointed by bufferra. The dimensions of the texture are loaded into dimx and dimy */
        retval = load_ppm("selfi.ppm", &bufferra, &dimx, &dimy);
        if (retval !=1)
            {
            printf("Ez dago texturaren fitxategia (selfi.ppm)\n");
            exit(-1);
            }

	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        denak = 1;//1
        lineak = 0;//1
        objektuak = 1;//1
        foptr = 0;
        sel_ptr = 0;
        fkptr = 0;
        selk_ptr = 0;
        aldaketa = 'r';
        ald_lokala = 1;
        kamara = 0;
        ikuspegia = 1;
        analisia = 0;
        bektoreak = 0;
        back_culling = 0;
        proiekzioa = 0;
        gouraud = 0;
        argia_aukeratu = 0;
        num_lights = 4;
        Ia.r = 200.0;
        Ia.g = 75.0;
        Ia.b = 255.0;  
        perspektiba_proiekzioa(&(mperspektiba[0]), -5.0, 5.0, -5.0, 5.0, 5.0, 500.0);

        objektua_sortu = 1;
        read_from_file("cam.obj");
        fkptr->mptr->m[11] = 150;

        selk_ptr->rgb.b = 0.0;
        selk_ptr->rgb.b = 0.0;
        selk_ptr->rgb.b = 255.0;

        selk_ptr->Ka.r = 0.0;
        selk_ptr->Ka.g = 0.3;
        selk_ptr->Ka.b = 1.0;
        
        selk_ptr->kd.r = 0.5922;
        selk_ptr->kd.g = 0.0166;
        selk_ptr->kd.b = 0.0;

        selk_ptr->ks.r = 0.5974;
        selk_ptr->ks.g = 0.2084;
        selk_ptr->ks.b = 0.2084;

        selk_ptr->ns = 100.2237;
        

        objektua_sortu = 0;
        
        read_from_file("x_wing.obj");
        sel_ptr->mptr->m[3] = -100;

        sel_ptr->rgb.r = 102.0; 
        sel_ptr->rgb.g = 204.0; 
        sel_ptr->rgb.b = 76.5; 

        sel_ptr->Ka.r = 0.4; 
        sel_ptr->Ka.g = 0.8;
        sel_ptr->Ka.b = 0.3;

        sel_ptr->kd.r = 0.3;
        sel_ptr->kd.g = 0.4;
        sel_ptr->kd.b = 0.6;

        sel_ptr->ks.r = 0.3;
        sel_ptr->ks.g = 0.1;
        sel_ptr->ks.b = 0.2;

        sel_ptr->ns = 60;
        
        read_from_file("r_falke.obj");
        sel_ptr->mptr->m[3] = 100;

        sel_ptr->rgb.r = 10.0;
        sel_ptr->rgb.g = 255.0;
        sel_ptr->rgb.b = 130.0;

        sel_ptr->Ka.r = 0.1; 
        sel_ptr->Ka.g = 1.0;
        sel_ptr->Ka.b = 0.2;

        sel_ptr->kd.r = 0.8;
        sel_ptr->kd.g = 0.8;
        sel_ptr->kd.b = 0.8;

        sel_ptr->ks.r = 0.5;
        sel_ptr->ks.g = 0.5;
        sel_ptr->ks.b = 0.5;

        sel_ptr->ns = 225.0;
        
        read_from_file("cam.obj");
        sel_ptr->rgb.r = 255.0;
        sel_ptr->rgb.g = 215.0;
        sel_ptr->rgb.b = 0.0;

        sel_ptr->Ka.r = 0.8;
        sel_ptr->Ka.g = 0.72;
        sel_ptr->Ka.b = 0.0;

        sel_ptr->kd.r = 0.8;
        sel_ptr->kd.g = 0.6;
        sel_ptr->kd.b = 0.4;

        sel_ptr->ks.r = 0.9;
        sel_ptr->ks.g = 0.8;
        sel_ptr->ks.b = 0.5;

        sel_ptr->ns = 80.0;
        
        if (argc>1) read_from_file(argv[1]);

        argiak_sortu();

	glutMainLoop();
	return 0;
}

