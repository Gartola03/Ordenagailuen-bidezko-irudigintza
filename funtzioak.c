#include <GL/glut.h>
#include <math.h>
#include "obj.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void print_matrizea(char *str, double *m){
    int i;

    printf("%s\n",str);
    for (i = 0;i<4;i++)
    printf("%lf, %lf, %lf, %lf\n", m[i*4], m[i*4+1], m[i*4+2], m[i*4+3]);
}

double modulua(double *v){
    return sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
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

double biderketa_eskalarra(double *v1, double *v2){
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void matrizen_eta_bektorearen_biderketa(double *v, double *m_ptr, double *k){
    int i;
    for (i=0; i<3; i++){
        v[i] = biderketa_eskalarra(&(m_ptr[i*4]), &(k[0]));
    }

}

void puntua_sortu(punto *p, vertex erpina){
    p->x = erpina.coord.x;
    p->y = erpina.coord.y;
    p->z = erpina.coord.z;
    p->u = erpina.u;
    p->v = erpina.v;
}

void biderketa_bektoriala(double *v, double *v1, double *v2){

    v[0] = v1[1]*v2[2] - v2[1]*v1[2];
    v[1] = -(v1[0]*v2[2] - v2[0]*v1[2]);
    v[2] = v1[0]*v2[1] - v2[0]*v1[1];

}

void bektorea_normalizatu(double *v){
    double mod;
    mod = modulua(&(v[0]));
    v[0] = v[0] / mod;
    v[1] = v[1] / mod;
    v[2] = v[2] / mod;
}

void bektore_normala_kalkulatu(double *nptr, punto p1, punto p2, punto p3){
    double v01[3];
    double v02[3];
    double v[3];

    v01[0] = p2.x-p1.x;
    v01[1] = p2.y-p1.y;
    v01[2] = p2.z-p1.z;

    v02[0] = p3.x-p1.x;
    v02[1] = p3.y-p1.y;
    v02[2] = p3.z-p1.z;

    biderketa_bektoriala(&(nptr[0]), &(v01[0]), &(v02[0]));
    bektorea_normalizatu(&(nptr[0]));

}

void bektore_normala_guztiak_kalkulatu(object3d *optr){
    int i, j, ind1, ind2, ind3, ind;
    double nlag[3];
    punto p1, p2, p3;
    
    for (i=0; i<optr->num_faces; i++){
        ind1 = optr->face_table[i].vertex_ind_table[0];
        ind2 = optr->face_table[i].vertex_ind_table[1];
        ind3 = optr->face_table[i].vertex_ind_table[2];
        
        puntua_sortu(&p1, optr->vertex_table[ind1]);
        puntua_sortu(&p2, optr->vertex_table[ind2]);
        puntua_sortu(&p3, optr->vertex_table[ind3]);

        bektore_normala_kalkulatu(&(optr->face_table[i].N[0]),p1,p2,p3);
        for (j=0; j <optr->face_table[i].num_vertices; j++){
            ind = optr->face_table[i].vertex_ind_table[j];
            optr->vertex_table[ind].N[0] += optr->face_table[i].N[0];
            optr->vertex_table[ind].N[1] += optr->face_table[i].N[1];
            optr->vertex_table[ind].N[2] += optr->face_table[i].N[2];
        }
    }
    
    for (i=0; i<optr->num_vertices; i++){
        bektorea_normalizatu(&(optr->vertex_table[i].N[0]));
    }
    
}

void erreferentzia_sistemaren_aldaketa(double *mesa, object3d *optr){
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

void modelview_kalkulatu(double *modelview, double *mesa, object3d *auxptr){
    double mlag[16];
    biderkatu_matrizeak(&(modelview[0]),&(mesa[0]),&(auxptr->mptr->m[0]));
}

void kamara_arreta_puntuari_begiratu(object3d *selk_ptr, object3d *sel_ptr){
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
    bektorea_normalizatu(&(z[0]));

    biderketa_bektoriala(&(x[0]), &(vup[0]), &(z[0]));
    bektorea_normalizatu(&(x[0]));

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

void biraketa_matrizea_analisi_moduan(object3d *selk_ptr, object3d *sel_ptr , double *mlag, char biraketa, double sinalfa, double cosalfa){
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

void perspektiba_proiekzioa(double *mperspektiba, double l, double r, double b, double t, double n, double f){
    int i;
    for(i=0; i<16; i++) mperspektiba[i] = 0.0;
    mperspektiba[0] = 2*n/(r-l);
    mperspektiba[2] = (r+l)/(r-l);
    mperspektiba[5] = 2*n/(t-b);;
    mperspektiba[6] = (t+b)/(t-b);
    mperspektiba[10] = (-f+n)/(f-n);
    mperspektiba[11] = -2*f*n/(f-n);
    mperspektiba[14] = -1;

}

void mxp(punto *pptr, double m[16], punto p){
    pptr->x = m[0]*p.x + m[1]*p.y + m[2]*p.z + m[3];
    pptr->y = m[4]*p.x + m[5]*p.y + m[6]*p.z + m[7];
    pptr->z = m[8]*p.x + m[9]*p.y + m[10]*p.z + m[11];
    pptr->u = p.u;
    pptr->v = p.v;

}

int mxpp(punto *pptr, double m[16], punto p){
    double x1,y1,z1,w;

    x1 = m[0]*p.x + m[1]*p.y + m[2]*p.z + m[3];
    y1 = m[4]*p.x + m[5]*p.y + m[6]*p.z + m[7];
    z1 = m[8]*p.x + m[9]*p.y + m[10]*p.z + m[11];
    w = m[12]*p.x + m[13]*p.y + m[14]*p.z + m[15];

    if (w <= 0){
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

void objektuari_aldaketa_sartu_ezk(double m[16], object3d *optr){
    double *mptrm;
    mlist *bmptr;

    mptrm=&(optr->mptr->m[0]);
    bmptr=(mlist*)malloc(sizeof(mlist));
    biderkatu_matrizeak(&(bmptr->m[0]),&(m[0]),mptrm);
    bmptr->hptr = &(optr->mptr[0]);
    optr->mptr = &(bmptr[0]);
}

void objektuari_aldaketa_sartu_esk(double m[16], object3d *optr){
    double *mptrm;
    mlist *bmptr;

    mptrm=&(optr->mptr->m[0]);
    bmptr=(mlist*)malloc(sizeof(mlist));
    biderkatu_matrizeak(&(bmptr->m[0]),mptrm, &(m[0]));
    bmptr->hptr = &(optr->mptr[0]);
    optr->mptr = &(bmptr[0]);
}

void gouraud_bektoreak_marraztu(object3d *optr, double *modelview, double *mperspektiba, int proiekzioa, int indexx){
    double n[3], pmodelview[16];
    punto p;
    int i,a;
    if (indexx != -1){
        puntua_sortu(&p, optr->vertex_table[indexx]);
        mxp(&p,modelview,p);

        if (proiekzioa == 0){
            a = mxpp(&p,mperspektiba,p);
            if (a == 1) return;
            matrizen_eta_bektorearen_biderketa(&(n[0]),&(mperspektiba[0]) ,&(optr->vertex_table[indexx].Ncam[0]));
            glBegin(GL_LINES);
                glVertex3d(p.x, p.y, p.z);
                glVertex3d(p.x + 50* n[0], p.y + 50* n[1], p.z + 50* n[2]);
            glEnd();

        } else{
            glBegin(GL_LINES);
                glVertex3d(p.x, p.y, p.z);
                glVertex3d(p.x + 50* optr->vertex_table[indexx].Ncam[0], p.y + 50* optr->vertex_table[indexx].Ncam[1], p.z + 50* optr->vertex_table[indexx].Ncam[2]);
            glEnd();
        }

    }else{
        for (i = 0; i<optr->num_vertices; i++){
            puntua_sortu(&p, optr->vertex_table[i]);
            mxp(&p,modelview,p);

            if (proiekzioa == 0){
                a = mxpp(&p,mperspektiba,p);
                if (a == 1) return;
                matrizen_eta_bektorearen_biderketa(&(n[0]),&(mperspektiba[0]) ,&(optr->vertex_table[i].Ncam[0]));
                glBegin(GL_LINES);
                    glVertex3d(p.x, p.y, p.z);
                    glVertex3d(p.x + 50* n[0], p.y + 50* n[1], p.z + 50* n[2]);
                glEnd();

            } else{
                glBegin(GL_LINES);
                    glVertex3d(p.x, p.y, p.z);
                    glVertex3d(p.x + 50* optr->vertex_table[i].Ncam[0], p.y + 50* optr->vertex_table[i].Ncam[1], p.z + 50* optr->vertex_table[i].Ncam[2]);
                glEnd();
            }   
        }
    }  
}

void flat_bektoreak_marraztu(object3d *optr, double *modelview, double *mperspektiba, int proiekzioa, int i){
    double n[3];
    punto p1;
    int ind1,a;

    ind1 = optr->face_table[i].vertex_ind_table[0];
    puntua_sortu(&p1, optr->vertex_table[ind1]);
    mxp(&p1,modelview,p1);

    if (proiekzioa == 0){
        a = mxpp(&p1,mperspektiba,p1);
        if (a == 1) return;
        matrizen_eta_bektorearen_biderketa(&(n[0]),&(mperspektiba[0]) ,&(optr->face_table[i].Ncam[0]));
        glBegin(GL_LINES);
            glVertex3d(p1.x, p1.y, p1.z);
            glVertex3d(p1.x + 50* n[0], p1.y + 50* n[1], p1.z + 50* n[2]);
        glEnd();

    } else{
        glBegin(GL_LINES);
            glVertex3d(p1.x, p1.y, p1.z);
            glVertex3d(p1.x + 50* optr->face_table[i].Ncam[0], p1.y + 50* optr->face_table[i].Ncam[1], p1.z + 50* optr->face_table[i].Ncam[2]);
        glEnd();
    }   
}

void light_bektore_cam(punto p1, light *light_table, double *modelview){
    //A1
    matrizen_eta_bektorearen_biderketa(&(light_table[0].camdir[0]), &(modelview[0]), &(light_table[0].dir[0]));
    bektorea_normalizatu(&(light_table[0].camdir[0]));
    //A2
    matrizen_eta_bektorearen_biderketa(&(light_table[1].camdir[0]), &(modelview[0]), &(light_table[1].pos[0]));
    light_table[1].camdir[0] -= p1.x;
    light_table[1].camdir[1] -= p1.y;
    light_table[1].camdir[2] -= p1.z; 
    bektorea_normalizatu(&(light_table[1].camdir[0]));

    //A3 = L = Nor(modelview -P1.cam)
    matrizen_eta_bektorearen_biderketa(&(light_table[2].camdir[0]), &(modelview[0]), &(light_table[2].pos[0]));
    light_table[2].camdir[0] -= p1.x;
    light_table[2].camdir[1] -= p1.y;
    light_table[2].camdir[2] -= p1.z; 
    bektorea_normalizatu(&(light_table[2].camdir[0]));
    //A4 = L = -P1.cam
    light_table[3].camdir[0] = -p1.x;
    light_table[3].camdir[1] = -p1.y;
    light_table[3].camdir[2] = -p1.z;
    bektorea_normalizatu(&(light_table[3].camdir[0]));
}

void h_kalkulatu(double *H, double *vkamara, light *light_table){
    double l[3], mod;
    int j, num_lights;

    num_lights = 4;
    for (j = 0; j<num_lights; j++){
        H[j*3] = vkamara[0] + light_table[j].camdir[0];
        H[j*3 + 1] = vkamara[1] + light_table[j].camdir[1];
        H[j*3 + 2] = vkamara[2] + light_table[j].camdir[2];
        mod = sqrt(pow(H[j*3], 2) + pow(H[j*3 + 1], 2) + pow(H[j*3 + 2], 2));
        H[j*3] = H[j*3]/mod;
        H[j*3 + 1] = H[j*3 + 1]/mod;
        H[j*3 + 2] = H[j*3 + 2]/mod;
    }
}
 
void argi_kalkulua(color3 *I, object3d *optr, color3 Ia, int aurp, int gouraud, int ind, light *light_table, double *H){
    int i;
    double nl, nh, nhns, nf;
    double h[3], vn[3];

    I->r = Ia.r * optr->Ka.r;
    I->g = Ia.g * optr->Ka.g;
    I->b = Ia.b * optr->Ka.b;

    if (gouraud == 0){
        for (i = 0; i<2; i++){
            if (light_table[i].onoff == 1){
                nl = biderketa_eskalarra(&(optr->face_table[aurp].Ncam[0]), &(light_table[i].camdir[0]));
                if (nl < 0) nl = 0;
                h[0] = H[3*i];
                h[1] = H[3*i+1];
                h[2] = H[3*i+2];
                nh = biderketa_eskalarra(&(optr->face_table[aurp].Ncam[0]), &(h[0]));
                if (nh < 0) nh = 0;
                I->r += light_table[i].I.r * (optr->kd.r * nl + optr->ks.r * pow(nh,optr->ns));
                I->g += light_table[i].I.g * (optr->kd.g * nl + optr->ks.g * pow(nh,optr->ns));
                I->b += light_table[i].I.b * (optr->kd.b * nl + optr->ks.b * pow(nh,optr->ns)); 
            }
        }
        for (i = 2; i<4; i++){
            vn[0] = -light_table[i].camdir[0];
            vn[1] = -light_table[i].camdir[1];
            vn[2] = -light_table[i].camdir[2];
            nf = biderketa_eskalarra(&(light_table[i].foko[0]), &(vn[0]));
            
            if (light_table[i].onoff == 1  && nf > cos(light_table[i].aperture)){
                nl = biderketa_eskalarra(&(optr->face_table[aurp].Ncam[0]), &(light_table[i].camdir[0]));
                if (nl < 0) nl = 0;
                h[0] = H[3*i];
                h[1] = H[3*i+1];
                h[2] = H[3*i+2];
                nh = biderketa_eskalarra(&(optr->face_table[aurp].Ncam[0]), &(h[0]));
                if (nh < 0) nh = 0;
                I->r += light_table[i].I.r * (optr->kd.r * nl + optr->ks.r * pow(nh,optr->ns));
                I->g += light_table[i].I.g * (optr->kd.g * nl + optr->ks.g * pow(nh,optr->ns));
                I->b += light_table[i].I.b * (optr->kd.b * nl + optr->ks.b * pow(nh,optr->ns)); 
            }
        }
    }else{
        for (i = 0; i<2; i++){
            if (light_table[i].onoff == 1){
                nl = biderketa_eskalarra(&(optr->vertex_table[ind].Ncam[0]), &(light_table[i].camdir[0]));
                if (nl < 0) nl = 0;
                h[0] = H[3*i];
                h[1] = H[3*i+1];
                h[2] = H[3*i+2];
                nh = biderketa_eskalarra(&(optr->vertex_table[ind].Ncam[0]), &(h[0]));
                if (nh < 0) nh = 0;
                I->r += light_table[i].I.r * (optr->kd.r * nl + optr->ks.r * pow(nh,optr->ns));
                I->g += light_table[i].I.g * (optr->kd.g * nl + optr->ks.g * pow(nh,optr->ns));
                I->b += light_table[i].I.b * (optr->kd.b * nl + optr->ks.b * pow(nh,optr->ns)); 
            }
        }
        for (i = 2; i<4; i++){
            vn[0] = -light_table[i].camdir[0];
            vn[1] = -light_table[i].camdir[1];
            vn[2] = -light_table[i].camdir[2];
            nf = biderketa_eskalarra(&(light_table[i].foko[0]), &(vn[0]));
        
            if (light_table[i].onoff == 1 && nf > cos(light_table[i].aperture)){
                nl = biderketa_eskalarra(&(optr->vertex_table[ind].Ncam[0]), &(light_table[i].camdir[0]));
                if (nl < 0) nl = 0;
                h[0] = H[3*i];
                h[1] = H[3*i+1];
                h[2] = H[3*i+2];
                nh = biderketa_eskalarra(&(optr->vertex_table[ind].Ncam[0]), &(h[0]));
                if (nh < 0) nh = 0;
                I->r += light_table[i].I.r * (optr->kd.r * nl + optr->ks.r * pow(nh,optr->ns));
                I->g += light_table[i].I.g * (optr->kd.g * nl + optr->ks.g * pow(nh,optr->ns));
                I->b += light_table[i].I.b * (optr->kd.b * nl + optr->ks.b * pow(nh,optr->ns)); 
            }
        }
    }
    if (I->r > 255.0) I->r = 255.0;
    if (I->g > 255.0) I->g = 255.0;
    if (I->b > 255.0) I->b = 255.0;
    
}