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
#include "cargar-triangulo.h"
#include <stdlib.h>
#include <math.h>

typedef struct mlist
    {
    double m[16]; 
    struct mlist *hptr;
    } mlist;
          
typedef struct triobj
    {
    hiruki *triptr;
    int num_triangles;
    unsigned char *kolorea;
    mlist *mptr;
    struct triobj *hptr;
    } triobj;

// testuraren informazioa
// información de textura

extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int * dimyptr);
unsigned char *bufferra;
int dimx,dimy;

int indexx;
hiruki *triangulosptr;
triobj *foptr;
triobj *sel_ptr;
triobj *fkptr;
triobj *selk_ptr;
int denak;
int lineak;
int objektuak;
char aldaketa;
int ald_lokala;
int kamara;
int kamara_sortu;
int ikuspegia;
int analisia;
int bektoreak;
int back_culling;
int proiekzioa;
double mesa[16];
double modelview1[16];
double mperspektiba[16];
double vkamara[3];

char fitxiz[100];

double modulua(double *v){
    return sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
}

double biderketa_bektoriala(double *v, double *v1, double *v2){

    v[0] = v1[1]*v2[2] - v2[1]*v1[2];
    v[1] = -(v1[0]*v2[2] - v2[0]*v1[2]);
    v[2] = v1[0]*v2[1] - v2[0]*v1[1];
    
}

void bektore_normala_kalkulatu(double *nptr, punto p1, punto p2, punto p3){
    double v01[3];
    double v02[3];
    double v[3];
    double mod;

    v01[0] = p2.x-p1.x;
    v01[1] = p2.y-p1.y;
    v01[2] = p2.z-p1.z;

    v02[0] = p3.x-p1.x;
    v02[1] = p3.y-p1.y;
    v02[2] = p3.z-p1.z;

    biderketa_bektoriala(&(v[0]), &(v01[0]), &(v02[0]));
    mod = modulua(&(v[0]));

    nptr[0] = v[0] / mod;
    nptr[1] = v[1] / mod;
    nptr[2] = v[2] / mod;

}
void perspektiba_proiekzioa(double l, double r, double b, double t, double n, double f){
    int i;
    for(i=0; i<16; i++) mperspektiba[i] = 0.0;
    mperspektiba[0] = 2*n/(r-l);
    mperspektiba[2] = (r+l)/(r-l);
    mperspektiba[5] = 2*n/(t-b);;
    mperspektiba[6] = (t+b)/(t-b);
    mperspektiba[10] = (-f+n)/(f-n);
    mperspektiba[11] = -2*f*n/(f-n);
    mperspektiba[14] =-1;

}

void kamara_arreta_puntuari_begiratu(){
    double mz, mx, my;
    double z[3];
    double x[3];
    double y[3];
    double k[3];
    double at[3];
    double vup[3];

    //Kamararen kokapena
    k[0] = selk_ptr->mptr->m[3];
    k[1] = selk_ptr->mptr->m[7];
    k[2] = selk_ptr->mptr->m[11];
    //Objekuaren arreta puntua
    at[0] = sel_ptr->mptr->m[3];
    at[1] = sel_ptr->mptr->m[7];
    at[2] = sel_ptr->mptr->m[11];
    //Vup bektorea
    vup[0] = selk_ptr->mptr->m[1]; 
    vup[1] = selk_ptr->mptr->m[5]; 
    vup[2] = selk_ptr->mptr->m[9]; 

    z[0] = k[0] - at[0];
    z[1] = k[1] - at[1];
    z[2] = k[2] - at[2];
    mz = modulua(&(z[0]));
    z[0] = z[0] / mz;
    z[1] = z[1] / mz;
    z[2] = z[2] / mz;

    biderketa_bektoriala(&(x[0]), &(vup[0]), &(z[0]));
    mx = modulua(&(x[0]));
    x[0] = x[0] / mx;
    x[1] = x[1] / mx;
    x[2] = x[2] / mx;

    biderketa_bektoriala(&(y[0]), &(z[0]), &(x[0]));

    selk_ptr->mptr->m[0] = x[0];
    selk_ptr->mptr->m[4] = x[1];
    selk_ptr->mptr->m[8] = x[2];

    selk_ptr->mptr->m[1] = y[0];
    selk_ptr->mptr->m[5] = y[1];
    selk_ptr->mptr->m[9] = y[2];

    selk_ptr->mptr->m[2] = z[0];
    selk_ptr->mptr->m[6] = z[1];
    selk_ptr->mptr->m[10] = z[2];

}

double biderketa_eskalarra(double *v1, double *v2){
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; 
}

void matrizen_eta_bektorearen_biderketa(double *v, double *m_ptr, double *k){
    int i;
    for (i=0; i<3; i++){
        v[i] = biderketa_eskalarra(&(m_ptr[i*4]), &(k[0]));
    }
    
}

void erreferentzia_sistemaren_aldaketa(triobj *optr){
    double *maptr; 
    double xc[3];
    double yc[3];
    double zc[3];
    double k[3];
    double ex, ey, ez;

    maptr = &(optr->mptr->m[0]);

    xc[0] = maptr[0];
    xc[1] = maptr[4];
    xc[2] = maptr[8];
    
    yc[0] = maptr[1];
    yc[1] = maptr[5];
    yc[2] = maptr[9];
    
    zc[0] = maptr[2];
    zc[1] = maptr[6];
    zc[2] = maptr[10];

    k[0] = maptr[3];
    k[1] = maptr[7];
    k[2] = maptr[11];

    ex = biderketa_eskalarra(&(k[0]), &(xc[0]));
    ey = biderketa_eskalarra(&(k[0]), &(yc[0]));
    ez = biderketa_eskalarra(&(k[0]), &(zc[0]));

    mesa[0] = xc[0];
    mesa[1] = xc[1];
    mesa[2] = xc[2];

    mesa[4] = yc[0];
    mesa[5] = yc[1];
    mesa[6] = yc[2];

    mesa[8] = zc[0];
    mesa[9] = zc[1];
    mesa[10] = zc[2];

    mesa[3] = -ex;
    mesa[7] = -ey;
    mesa[11] = -ez;

    mesa[12] = 0.0;
    mesa[13] = 0.0;
    mesa[14] = 0.0;

    mesa[15] = 1.0;

}

void biderkatu_matrizeak(double *emaitzaptr, double *eskerptr, double *eskuinptr){
    int i, j, k;
    double emaitza;
    for (i=0; i<4; i++){
        for(j=0; j<4; j++){
            emaitza=0.0;
            for(k=0; k<4; k++){
                emaitza+= eskerptr[4*i+k]*eskuinptr[4*k+j];
            }
            emaitzaptr[4*i+j]=emaitza;
        }
    }
}

void biraketa_matrizea_analisi_moduan(double *mlag, char biraketa, double sinalfa, double cosalfa){
    double atx, aty, atz, x, y, z;
    double mbira[16];
    double mr[16];
    double mtat[16];
    double mt_at[16];
    int i;

    for (i=0; i<16; i++){
        mtat[i] = 0.0;
        mt_at[i] = 0.0;
        mlag[i] = 0.0;
    } 
    mtat[0] = 1.0;
    mtat[5] = 1.0;
    mtat[10] = 1.0;
    mtat[15] = 1.0;

    mt_at[0] = 1.0;
    mt_at[5] = 1.0;
    mt_at[10] = 1.0;
    mt_at[15] = 1.0;

    mlag[0] = 1.0;
    mlag[5] = 1.0;
    mlag[10] = 1.0;
    mlag[15] = 1.0;

    
    atx = sel_ptr->mptr->m[3];
    aty = sel_ptr->mptr->m[7];
    atz = sel_ptr->mptr->m[11];
    
    mtat[3] = atx;
    mtat[7] = aty;
    mtat[11] = atz;

    mt_at[3] = -atx;
    mt_at[7] = -aty;
    mt_at[11] = -atz;

    switch(biraketa){
        case 'x':
            x = selk_ptr->mptr->m[0]; 
            y = selk_ptr->mptr->m[4]; 
            z = selk_ptr->mptr->m[8]; 
            break;
        case 'y':
            x = selk_ptr->mptr->m[1]; 
            y = selk_ptr->mptr->m[5];
            z = selk_ptr->mptr->m[9];
            break;
        case 'z':
            x = selk_ptr->mptr->m[2]; 
            y = selk_ptr->mptr->m[6]; 
            z = selk_ptr->mptr->m[10]; 
            break;
    }

    mr[0] = cosalfa + (1-cosalfa)*pow(x,2);
    mr[1] = (1-cosalfa)*x*y - z*sinalfa;
    mr[2] = (1-cosalfa)*x*z + y*sinalfa;
    mr[3] = 0.0;

    mr[4] = (1-cosalfa)*x*y + z*sinalfa;
    mr[5] = cosalfa + (1-cosalfa)*pow(y,2);
    mr[6] = (1-cosalfa)*y*z - x*sinalfa;
    mr[7] = 0.0;

    mr[8] = (1-cosalfa)*x*z - y*sinalfa;
    mr[9] = (1-cosalfa)*y*z + x*sinalfa;
    mr[10] = cosalfa + (1-cosalfa)*pow(z,2);
    mr[11] = 0.0;

    mr[12] = 0.0;
    mr[13] = 0.0;
    mr[14] = 0.0;
    mr[15] = 1.0;

    biderkatu_matrizeak(&(mbira[0]),&(mr[0]), &(mt_at[0]));
    biderkatu_matrizeak(&(mlag[0]),&(mtat[0]), &(mbira[0]));

}


void objektuari_aldaketa_sartu_ezk(double m[16], triobj *optr)
{
    double *mptrm;
    mlist *bmptr;
    
    mptrm=&(optr->mptr->m[0]);
    bmptr=(mlist*)malloc(sizeof(mlist));
    biderkatu_matrizeak(&(bmptr->m[0]),&(m[0]),mptrm);
    bmptr->hptr = &(optr->mptr[0]);
    optr->mptr = &(bmptr[0]);
}



void objektuari_aldaketa_sartu_esk(double m[16], triobj *optr)
{
    double *mptrm;
    mlist *bmptr;
    
    mptrm=&(optr->mptr->m[0]);
    bmptr=(mlist*)malloc(sizeof(mlist));
    biderkatu_matrizeak(&(bmptr->m[0]),mptrm, &(m[0]));
    bmptr->hptr = &(optr->mptr[0]);
    optr->mptr = &(bmptr[0]);
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
void  dibujar_linea_z(int linea,float c1x, float c1z, float c1u,float c1v,float c2x,float c2z,float c2u,float c2v){
    float xkoord,zkoord;
    float u,v;
    float t,deltat;
    unsigned char r,g,b;
    unsigned char *colorv;

    glBegin( GL_POINTS );
    if(c1x != c2x){
        deltat = 1.0/(c2x-c1x);
        //Lerroan xkoord, zkoord, u eta v posizio berriak eguneratzen du. Zuzenaren koordenatu barizentrikoa erabiltzen du
        for (t = 1.0; t>=0;t-=deltat)
        {
            xkoord = c1x*t +(1-t)*c2x;
            zkoord = c1z*t +(1-t)*c2z;
            u = c1u*t +(1-t)*c2u;
            v = c1v*t +(1-t)*c2v;
            //Texturaren posizioa jasotzen du
            colorv=  color_textura(u, v); 
            r= colorv[0];
            g=colorv[1];
            b=colorv[2];    
            glColor3ub(r,g,b);
            glVertex3f(xkoord, linea, zkoord );
        }
    }

    glEnd();
}


void print_matrizea(char *str, triobj *auxptr){
    int i;

    printf("%s\n",str);
    for (i = 0;i<4;i++)
    printf("%lf, %lf, %lf, %lf\n", auxptr->mptr->m[i*4], auxptr->mptr->m[i*4+1], auxptr->mptr->m[i*4+2], auxptr->mptr->m[i*4+3]);
}


void info (char *str){
    printf("\n");
    printf("INFORMAZIOA:\n");
    if (proiekzioa == 1){
        printf("Proiekzioa paraleloan ikusten ari zara.\n");
    } else{
        printf("Proiekzioa perspektiban ikusten ari zara.\n");
    }

    if (ikuspegia == 1){
        printf("Kamara ikuspegian ikusten ari zara.\n");
    } else{
        printf("Objektuaren ikuspegian ikusten ari zara.\n");
    }

    if (kamara == 1){
        printf("Kamara aukeratzen ari zara.\n");
    } else{
        printf("Objektua aukeratzen ari zara.\n");
    }

    if (ald_lokala == 1){
        printf("Aldaketa lokala (hegan moduan) eragiten ari zara.\n");
    } else{
        printf("Aldaketa globala (analisi moduan) eragiten ari zara.\n");
    }

    if (aldaketa == 't'){
        printf("Translaziok eragingo dizkiozu.\n");
    } else{
        printf("Biraketak eragingo dizkiozu.\n");
    }

    if (bektoreak == 1){
        printf("Bektoreak aktibatuta daude.\n");
    } else{
        printf("Bektoreak desakibatuta daude.\n");
    }

    if (back_culling == 0){
        printf("Atze-aurpegien marrazketa egiten ari zara.\n");
    } else{
        printf("Atze-aurpegien ezabaketa egiten ari zara.\n");
    }
    printf("\n");
}

int mxp(punto *pptr, double m[16], punto p)
{
    pptr->x = m[0]*p.x + m[1]*p.y + m[2]*p.z + m[3];
    pptr->y = m[4]*p.x + m[5]*p.y + m[6]*p.z + m[7];
    pptr->z = m[8]*p.x + m[9]*p.y + m[10]*p.z + m[11];
    pptr->u = p.u;
    pptr->v = p.v;

}

int mxpp(punto *pptr, double m[16], punto p)
{
    double x1,y1,z1,w;
    x1 = m[0]*p.x + m[1]*p.y + m[2]*p.z + m[3];
    y1 = m[4]*p.x + m[5]*p.y + m[6]*p.z + m[7];
    z1 = m[8]*p.x + m[9]*p.y + m[10]*p.z + m[11];
    w = m[12]*p.x + m[13]*p.y + m[14]*p.z + m[15];
    
    if (w == 0 || w ==-0){
        return 1;
    }
    x1 = 500*x1/w;
    y1 = 500*y1/w;
    z1 = -500*z1/w;

    pptr->x = x1;
    pptr->y = y1;
    pptr->z = z1;
    pptr->u = p.u;
    pptr->v = p.v;

    return 0;
}

int ikusten_da(triobj *aukptr, triobj *optr, int i){
    punto p1,p2,p3;
    hiruki *tptr;
    int k;

    tptr = optr->triptr+i;
    if (proiekzioa == 1){
        vkamara[0] = aukptr->mptr->m[2];
        vkamara[1] = aukptr->mptr->m[6];
        vkamara[2] = aukptr->mptr->m[10];
 
    } else{
        vkamara[0] = aukptr->mptr->m[3];
        vkamara[1] = aukptr->mptr->m[7];
        vkamara[2] = aukptr->mptr->m[11];

        mxp(&p1,optr->mptr->m,tptr->p1);

        vkamara[0] -= p1.x;
        vkamara[1] -= p1.y;
        vkamara[2] -= p1.z;
    }
    
    mxp(&p1,optr->mptr->m,tptr->p1);
    mxp(&p2,optr->mptr->m,tptr->p2);
    mxp(&p3,optr->mptr->m,tptr->p3);


    bektore_normala_kalkulatu(&(tptr->N[0]), p1, p2, p3);
    //marraztu edo ez
    glColor3ub(255, 255, 255);
    if (biderketa_eskalarra(&(tptr->N[0]), &(vkamara[0])) <= 0){
        if (back_culling == 1){
            return 1;
        }else{
            glColor3ub(255, 0, 0);
        }
    } 
    return 0;
}

void modelview_kalkulatu(triobj *auxptr){
    double mlag[16];
    //Modelview kalkulatzeko Modelview = Mesa * Matrizea
    if (proiekzioa == 1){
        biderkatu_matrizeak(&(modelview1[0]),&(mesa[0]),&(auxptr->mptr->m[0]));
    } else{
        biderkatu_matrizeak(&(mlag[0]),&(mesa[0]),&(auxptr->mptr->m[0]));
        biderkatu_matrizeak(&(modelview1[0]),&(mperspektiba[0]),&(mlag[0]));
    }
}

//Triangeluan lerroan marrazteko behar duen bi ebakidura puntuak kalkulatzen du. 
void ebaketa_kalkulatu(punto *Gptr, punto *Bptr, int h, punto *eptr)
{
    float alpha, beta;
    if(Gptr->y != Bptr ->y){
        alpha = (float)(Gptr->y-h)/ (Gptr->y - Bptr->y);
        beta = 1.0-alpha;
        eptr -> y = h;
        eptr -> x = beta * Gptr->x + alpha * Bptr->x;
        eptr -> z = beta * Gptr->z + alpha * Bptr->z;
        eptr -> u = beta * Gptr->u + alpha * Bptr->u;
        eptr -> v = beta * Gptr->v + alpha * Bptr->v;
        
    } 
}

void dibujar_triangulo(triobj *optr, int i){
    hiruki *tptr;
    triobj *auxptr;
    punto *pgoiptr, *pbeheptr, *perdiptr, *pezkeptr, *perdikoptr, *peskuinptr;
    float x1,h1,z1,u1,v1,x2,h2,z2,u2,v2,x3,h3,z3,u3,v3;
    float c1x,c1z,c1u,c1v,c2x,c2z,c2u,c2v;
    int linea;
    float cambio1,cambio1z,cambio1u,cambio1v,cambio2,cambio2z,cambio2u,cambio2v;
    punto p1,p2,p3;
    punto e1,e2,e;

    if (i >= optr->num_triangles) return;
    tptr = optr->triptr+i;


    if (proiekzioa == 1){
        mxp(&p1,modelview1,tptr->p1);
        mxp(&p2,modelview1,tptr->p2);
        mxp(&p3,modelview1,tptr->p3);

    } else{
        if (mxpp(&p1,modelview1,tptr->p1) == 1) return;
        if (mxpp(&p2,modelview1,tptr->p2) == 1) return;
        if (mxpp(&p3,modelview1,tptr->p3) == 1) return;
    }

    bektore_normala_kalkulatu(&(tptr->N[0]), p1, p2, p3);
    if (lineak == 1){
        glBegin(GL_POLYGON);
            glVertex3d(p1.x, p1.y, p1.z);
            glVertex3d(p2.x, p2.y, p2.z);
            glVertex3d(p3.x, p3.y, p3.z);
        glEnd();
        if (bektoreak == 1){
            glBegin(GL_LINES);
                glVertex3d(p1.x, p1.y, p1.z);
                glVertex3d(p1.x + 50* tptr->N[0], p1.y + 50* tptr->N[1], p1.z + 50* tptr->N[2]);
            glEnd();

        }

        return;
    }
    //  else 
    // TODO 
    // hemen azpikoa kendu eta triangelua testurarekin marrazten duen kodea sartu.
    // lo que sigue aqui hay que sustituir por el código adecuado que dibuja el triangulo con textura   

    //Triangeluaren 3 puntuen edo eskinen posizioak goiko, erdiko eta beheko puntuen arabera sailkatzen ditu
    if (p1.y > p2.y){
        pgoiptr = &p1;
        pbeheptr = &p2;
    } else{
        pgoiptr = &p2;
        pbeheptr = &p1;
    }

    if (p3.y > pgoiptr -> y){
        perdiptr = pgoiptr;
        pgoiptr = &p3;
    }else{
        if(p3.y < pbeheptr -> y){
            perdiptr = pbeheptr;
            pbeheptr = &p3;
        } 
        else perdiptr = &p3;
    }



    //Triangeluaren goiko eta erdikoaren puntuen artean, lerroak marrazten ditu
    if (pgoiptr->y != perdiptr->y)
    {
        for (i=pgoiptr -> y-1; i>perdiptr -> y; i--)
        {
            ebaketa_kalkulatu(pgoiptr, perdiptr, i, &e1);
            ebaketa_kalkulatu(pgoiptr, pbeheptr, i, &e2); 
            if (e1.x >e2.x){
                dibujar_linea_z(i, e2.x, e2.z, e2.u, e2.v, e1.x, e1.z, e1.u, e1.v);
            } else{
                dibujar_linea_z(i, e1.x, e1.z, e1.u, e1.v, e2.x, e2.z, e2.u, e2.v);
            } 

        
        }
    }

    //Triangeluaren erdiko eta behekoaren puntuen artean, lerroak marrazten ditu
    if (perdiptr->y != pbeheptr->y){
        for (i=perdiptr -> y-1; i>pbeheptr -> y; i--)
        {
            ebaketa_kalkulatu(perdiptr, pbeheptr, i, &e1);
            ebaketa_kalkulatu(pgoiptr, pbeheptr, i, &e2); 
            if (e1.x >e2.x){
                dibujar_linea_z(i,e2.x, e2.z, e2.u, e2.v, e1.x, e1.z, e1.u, e1.v);
            } else{
                dibujar_linea_z(i,e1.x, e1.z, e1.u, e1.v, e2.x, e2.z, e2.u, e2.v);
            } 
        }
    }
    if(pgoiptr->y == pbeheptr->y && pgoiptr->y == perdiptr->y){
        if(p1.x > p2.x){
            perdikoptr = &p1;
            pezkeptr = &p2;
        } else{
            perdikoptr = &p2;
            pezkeptr = &p1;
        }

        if(p3.x > perdikoptr -> x){
            peskuinptr = &p3;
        }
        else{
            if(p3.x < pezkeptr -> x){
                perdikoptr = pezkeptr;
                pezkeptr = &p3;
            }
        }
        dibujar_linea_z(perdikoptr -> y, pezkeptr->x, pezkeptr->z, pezkeptr->u, pezkeptr->v, peskuinptr->x, peskuinptr->z, peskuinptr->u, peskuinptr->v );
    } 
}




static void marraztu(void){
    float u,v;
    int i,j;
    triobj *auxptr, *aukptr;

    /*
    unsigned char* colorv;
    unsigned char r,g,b;
    */

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

    info("");
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

    erreferentzia_sistemaren_aldaketa(aukptr);

    if (fkptr != 0){       
        triangulosptr = selk_ptr->triptr;
        if (objektuak == 1){
            if (denak == 1){
                for (auxptr =fkptr; auxptr != 0; auxptr = auxptr->hptr){
                    modelview_kalkulatu(auxptr);
                    for (i =0; i < auxptr->num_triangles; i++){
                        if (ikusten_da(aukptr, auxptr, i) == 0) dibujar_triangulo(auxptr,i);                       
                    }
                }     
            } else{
                modelview_kalkulatu(selk_ptr);
                for (i =0; i < selk_ptr->num_triangles; i++)
                    {
                        if (ikusten_da(aukptr, selk_ptr, i) == 0) dibujar_triangulo(selk_ptr,i);
                    }
                }
        } else{
            modelview_kalkulatu(selk_ptr);
            if (ikusten_da(aukptr, selk_ptr, i) == 0){
                dibujar_triangulo(selk_ptr,indexx);
            }
        }
    } 

    // marrazteko objektuak behar dira
    // no se puede dibujar sin objetos
    if (foptr != 0){       
        triangulosptr = sel_ptr->triptr;
        if (objektuak == 1){
            if (denak == 1){
                for (auxptr =foptr; auxptr != 0; auxptr = auxptr->hptr){
                    modelview_kalkulatu(auxptr);
                    for (i =0; i < auxptr->num_triangles; i++){
                        if (ikusten_da(aukptr, auxptr, i) == 0){
                            dibujar_triangulo(auxptr,i);
                        }                          
                    }
                }
            } else{
                modelview_kalkulatu(sel_ptr);
                for (i =0; i < sel_ptr->num_triangles; i++){
                    if (ikusten_da(aukptr, sel_ptr, i) == 0){
                        dibujar_triangulo(sel_ptr,i);
                    }
                }
            }
        } else{
            modelview_kalkulatu(sel_ptr);
            if (ikusten_da(aukptr, sel_ptr, i) == 0){
                dibujar_triangulo(sel_ptr,indexx);
            }
        }
    } 

    glFlush();
}


void read_from_file(char *fitx){
    int i,n,retval;
    triobj *optr;
    hiruki *tptr;

    //printf("%s fitxategitik datuak hartzera\n",fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triptr), &(optr->kolorea));
    if (retval !=1 && retval != 15 && retval != 9) {
        printf("%s fitxategitik datuak hartzerakoan arazoak izan ditut\n    Problemas al leer\n",fitxiz);
        free(optr);
    } else {
        triangulosptr = optr->triptr;
        //printf("objektuaren matrizea...\n");
        optr->mptr = (mlist *)malloc(sizeof(mlist));
        for (i=0; i<16; i++) optr->mptr->m[i] =0.0;
        optr->mptr->m[0] = 1.0;
        optr->mptr->m[5] = 1.0;
        optr->mptr->m[10] = 1.0;
        optr->mptr->m[15] = 1.0;
        optr->mptr->hptr = 0;

        //Objektuaren triangelu guztien bektore normala kalkulatzen du
        for (n=0; n< optr->num_triangles; n++){
            tptr = optr->triptr+n;
            bektore_normala_kalkulatu(&(tptr->N[0]), tptr->p1, tptr->p2, tptr->p3);
        }
        
        //printf("objektu zerrendara doa informazioa...\n");
        if (kamara_sortu == 1){
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

void x_aldaketa(int dir){
    
    double p, alfa, sinalfa, cosalfa;
    triobj *auxptr;
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
        sinalfa = sin(-alfa);
        cosalfa = cos(-alfa);
    }
    if (kamara == 1){
        auxptr = selk_ptr;
        mlag[0]=cosalfa; 
        mlag[2]=sinalfa;
        mlag[8]=-sinalfa;
        mlag[10]=cosalfa;
        if (analisia == 1){
            biraketa = 'x';
            biraketa_matrizea_analisi_moduan(&(mlag[0]), biraketa, sinalfa, cosalfa);
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
            /*double mlag[16]={ 1,0,0,0
                            ,0,cos(alfa),-sin(alfa),0
                            ,0,sin(alfa),cos(alfa),0
                            ,0,0,0,1.0}; */
        }
    } 
    if(ald_lokala==1){
        objektuari_aldaketa_sartu_esk(mlag, auxptr);
    }else{
        objektuari_aldaketa_sartu_ezk(mlag, auxptr);
    } 
    
}


void y_aldaketa(int dir){
    double p, alfa, sinalfa, cosalfa;
    triobj *auxptr;
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
        sinalfa = sin(-alfa);
        cosalfa = cos(-alfa);
    }
    if (kamara == 1){
        auxptr = selk_ptr;
        mlag[5]=cosalfa;
        mlag[6]=-sinalfa;
        mlag[9]=sinalfa;
        mlag[10]=cosalfa;
        if (analisia == 1){
            biraketa = 'y';
            biraketa_matrizea_analisi_moduan(&(mlag[0]), biraketa, sinalfa, cosalfa);
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
            /*double mlag[16]={ cos(alfa),0,sin(alfa),0
                            ,0,1.0,0,0
                            ,-sin(alfa),0,cos(alfa),0
                            ,0,0,0,1.0}; */
        }
    }

    if(ald_lokala==1){
        objektuari_aldaketa_sartu_esk(mlag, auxptr);
    }else{
        objektuari_aldaketa_sartu_ezk(mlag, auxptr);
    }
     
    
}

void z_aldaketa(int dir){
    double p, alfa, sinalfa, cosalfa;
    triobj *auxptr;
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
        sinalfa = sin(-alfa);
        cosalfa = cos(-alfa);
    }
    mlag[11] = p; 
    if(aldaketa=='r'){
        mlag[0]=cosalfa; 
        mlag[1]=-sinalfa;
        mlag[4]=sinalfa;
        mlag[5]=cosalfa;
        mlag[11]=0.0; 
        /*double mlag[16]={ cos(alfa),-sin(alfa),0,0
                            ,sin(alfa),cos(alfa),0,0
                            ,0,0,1.0,0
                            ,0,0,0,1.0}; */
    } 
    
    if (kamara == 1 && analisia == 1){
        biraketa = 'z';
        biraketa_matrizea_analisi_moduan(&(mlag[0]), biraketa, sinalfa, cosalfa);
            
    }
        
    if(ald_lokala==1){
        objektuari_aldaketa_sartu_esk(mlag, auxptr);
    } else{
        objektuari_aldaketa_sartu_ezk(mlag, auxptr);
    }

}

void eskala_aldaketa(int dir){
    //proiekzioan, f aldatu
    double p;
    triobj *auxptr;
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
    triobj *auxptr;
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
    FILE *obj_file;

    switch(key){
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
                    kamara_arreta_puntuari_begiratu();
                }
            }
            else {
                ald_lokala = 1;
                analisia = 0; 
            }
            break;
        case 'c':
            if (kamara == 1) {
                kamara = 0; 
            }else {
                kamara = 1;
                if (analisia == 1){
                    kamara_arreta_puntuari_begiratu();
                }
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

        case 'f':
                /*Ask for file*/
                printf("idatzi fitxategi izena\n"); 
                scanf("%s", &(fitxiz[0]));
                if(kamara == 1){
                    kamara_sortu = 1;
                }else{
                    kamara_sortu = 0;
                }
                read_from_file(fitxiz);
                indexx = 0;
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
                if (kamara == 1){
                    selk_ptr = selk_ptr->hptr;
                    /*The selection is circular, thus if we move out of the list we go back to the first element*/
                    if (selk_ptr == 0) selk_ptr = fkptr;
                    indexx =0; // the selected polygon is the first one

                }else{
                    sel_ptr = sel_ptr->hptr;
                    /*The selection is circular, thus if we move out of the list we go back to the first element*/
                    if (sel_ptr == 0) sel_ptr = foptr;
                    indexx =0; // the selected polygon is the first one
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
            printf("Ez dago texturaren fitxategia (testura.ppm)\n");
            exit(-1);
            }
        
	glClearColor( 0.0f, 0.0f, 0.7f, 1.0f );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
        denak = 1;//1
        lineak =1;//1
        objektuak = 1;//1
        foptr = 0;
        sel_ptr = 0;
        aldaketa = 'r';
        ald_lokala = 1;
        fkptr = 0;
        selk_ptr = 0;
        kamara = 0;
        ikuspegia = 1;
        analisia = 0;
        bektoreak = 0;
        back_culling = 0;
        proiekzioa = 1; 
        perspektiba_proiekzioa(-5.0, 5.0, -5.0, 5.0, 5.0, 500.0);
        kamara_sortu = 1;
        read_from_file("cam.txt");
        fkptr -> mptr->m[11] = 200.0;
        kamara_sortu = 0;
        read_from_file("z.txt");
        foptr->mptr->m[3] = 200;
        kamara_sortu = 0;
        read_from_file("z.txt");
        foptr->mptr->m[3] = -200;
        if (argc>1) read_from_file(argv[1]);
	glutMainLoop();

	return 0;   
}
