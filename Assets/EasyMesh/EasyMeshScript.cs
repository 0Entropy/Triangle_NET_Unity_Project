using TriangleNet.Geometry;
using TriangleNet.Tools;
using UnityEngine;
using System;

public class EasyMeshScript : MonoBehaviour {

    int ON = 0;
    int OFF = -1;
    int WAIT = -2;
    int A = 3;
    int D = 4;
    int W = 5;

    int CLOSED = 0;
    int OPEN = 1;
    int INSIDE = 2;

    float PI = Mathf.PI;

    float SMALL = float.MinValue;
    float GREAT = float.MaxValue;

    const int MAX_NODES = 3000;

    char[] _name = new char[80];
    int len = 0;

    struct ele
    {
        public int i, j, k;
        public int ei, ej, ek;
        public int si, sj, sk;

        public int mark;             /* is it off (ON or OFF) */
        public int state;            /* is it (D)one, (A)ctive or (W)aiting */
        public int material;

        public float xv, yv, xin, yin, R, r, Det;

        public int new_numb;         /* used for renumeration */
    }

    ele[] elem = new ele[MAX_NODES * 2];

    struct sid
    {
        public int ea, eb;           /* left and right element */
        public int a, b, c, d;       /* left, right, start and end point */

        public int mark;             /* is it off, is on the boundary */

        public float s;

        public int new_numb;         /* used for renumeration */
    }
    sid[] side = new sid[MAX_NODES * 3];

    struct nod
    {
        public float x, y, F;

        public float sumx, sumy;
        public int Nne;

        public int mark;             /* is it off */

        public int next;             /* next node in the boundary chain */
        public int chain;            /* on which chains is the node */
        public int inserted;

        public int new_numb;         /* used for renumeration */
    }
    nod[] node = new nod[MAX_NODES];
    nod[] point = new nod[MAX_NODES / 2];

    struct seg
    {
        public int n0, n1;
        public int N;
        public int chain;
        public int bound;
        public int mark;
    }

    seg[] segment;// = new seg[];
 
    struct chai
    {
        public int s0, s1, type;
    }
    chai[] chain;// = new chai[];


    int Ne, Nn, Ns, Nc;             /* number of: elements, nodes, sides */
    int ugly;                       /* mora li biti globalna ??? */

    double area(nod na, nod nb, nod nc)
    {
     return 0.5d * (   ((nb).x-(na).x)*((nc).y-(na).y) 
		    - ((nb).y-(na).y)*((nc).x-(na).x));
    }

    double dist(nod na, nod nb)
    {
     return Mathf.Sqrt(((nb).x-(na).x)*((nb).x-(na).x) + ((nb).y-(na).y)*((nb).y-(na).y));
    }

    int in_elem(nod n)
    {
     int e;
 
     for(e=0; e<Ne; e++)    /* This must search through all elements ?? */
      {
       if(    area(n, node[elem[e].i], node[elem[e].j]) >= 0.0
           && area(n, node[elem[e].j], node[elem[e].k]) >= 0.0
           && area(n, node[elem[e].k], node[elem[e].i]) >= 0.0 )
   
       break;
      }
     return e;
    }

    void bowyer(int n, int spac)
    {
        int e,  s, swap; //i,
        nod vor = new nod();

         do  
          { 
           swap=0;
           for(s=0; s<Ns; s++)
           if(side[s].mark==0)
        /* if( !( (node[side[s].c].inserted>1 && node[side[s].d].bound==OFF && side[s].s<(node[side[s].c].F+node[side[s].d].F) ) ||
	          (node[side[s].d].inserted>1 && node[side[s].c].bound==OFF && side[s].s<(node[side[s].c].F+node[side[s].d].F) ) ) ) */
            {
             if(side[s].a==n)
              {e=side[s].eb; 
               if(e!=OFF)
	        {vor.x=elem[e].xv; 
	         vor.y=elem[e].yv;
	         if( dist(vor, node[n]) < elem[e].R )
	          {swap_side(s); swap=1;}}}
   
             else if(side[s].b==n)
              {e=side[s].ea; 
               if(e!=OFF)
	        {vor.x=elem[e].xv; 
	         vor.y=elem[e].yv;
	         if( dist(vor, node[n]) < elem[e].R )
	          {swap_side(s); swap=1;}}}
            }
          }
         while(swap==1);

       }

    /*=========================================================================*/
    void swap_side(int s)
    {
        int a, b, c, d, ea, eb;
        int eac = 0;
        int ead = 0;
        int ebc = 0;
        int ebd = 0;
        int sad = 0;
        int sac = 0;
        int sbc = 0;
        int sbd = 0;
        float sx, sy;

        ea = side[s].ea;
        eb = side[s].eb;
        a = side[s].a; b = side[s].b; c = side[s].c; d = side[s].d;

        if (elem[ea].ei == eb)
        {
            ead = elem[ea].ej; eac = elem[ea].ek;
            sad = elem[ea].sj; sac = elem[ea].sk;
        }

        if (elem[ea].ej == eb)
        {
            ead = elem[ea].ek; eac = elem[ea].ei;
            sad = elem[ea].sk; sac = elem[ea].si;
        }

        if (elem[ea].ek == eb)
        {
            ead = elem[ea].ei; eac = elem[ea].ej;
            sad = elem[ea].si; sac = elem[ea].sj;
        }

        if (elem[eb].ei == ea)
        {
            ebc = elem[eb].ej; ebd = elem[eb].ek;
            sbc = elem[eb].sj; sbd = elem[eb].sk;
        }

        if (elem[eb].ej == ea)
        {
            ebc = elem[eb].ek; ebd = elem[eb].ei;
            sbc = elem[eb].sk; sbd = elem[eb].si;
        }

        if (elem[eb].ek == ea)
        {
            ebc = elem[eb].ei; ebd = elem[eb].ej;
            sbc = elem[eb].si; sbd = elem[eb].sj;
        }

        elem[ea].i = a; elem[ea].j = b; elem[ea].k = d;
        elem[ea].ei = ebd; elem[ea].ej = ead; elem[ea].ek = eb;
        elem[ea].si = sbd; elem[ea].sj = sad; elem[ea].sk = s;

        elem[eb].i = a; elem[eb].j = c; elem[eb].k = b;
        elem[eb].ei = ebc; elem[eb].ej = ea; elem[eb].ek = eac;
        elem[eb].si = sbc; elem[eb].sj = s; elem[eb].sk = sac;

        if (eac != -1)
        {
            if (elem[eac].ei == ea) elem[eac].ei = eb;
            if (elem[eac].ej == ea) elem[eac].ej = eb;
            if (elem[eac].ek == ea) elem[eac].ek = eb;
        }

        if (ebd != -1)
        {
            if (elem[ebd].ei == eb) elem[ebd].ei = ea;
            if (elem[ebd].ej == eb) elem[ebd].ej = ea;
            if (elem[ebd].ek == eb) elem[ebd].ek = ea;
        }

        if (side[sad].ea == ea) { side[sad].a = b; }
        if (side[sad].eb == ea) { side[sad].b = b; }

        if (side[sbc].ea == eb) { side[sbc].a = a; }
        if (side[sbc].eb == eb) { side[sbc].b = a; }

        if (side[sbd].ea == eb) { side[sbd].ea = ea; side[sbd].a = a; }
        if (side[sbd].eb == eb) { side[sbd].eb = ea; side[sbd].b = a; }

        if (a < b)
        {
            side[s].c = a; side[s].d = b; side[s].a = d; side[s].b = c;
            side[s].ea = ea; side[s].eb = eb;
        }
        else
        {
            side[s].c = b; side[s].d = a; side[s].a = c; side[s].b = d;
            side[s].ea = eb; side[s].eb = ea;
        }

        sx = node[side[s].c].x - node[side[s].d].x;
        sy = node[side[s].c].y - node[side[s].d].y;
        side[s].s = Mathf.Sqrt(sx * sx + sy * sy);

        if (side[sac].ea == ea) { side[sac].ea = eb; side[sac].a = b; }
        if (side[sac].eb == ea) { side[sac].eb = eb; side[sac].b = b; }

        if (side[sad].ea == ea) { side[sad].a = b; }
        if (side[sad].eb == ea) { side[sad].b = b; }

        if (side[sbc].ea == eb) { side[sbc].a = a; }
        if (side[sbc].eb == eb) { side[sbc].b = a; }

        if (side[sbd].ea == eb) { side[sbd].ea = ea; side[sbd].a = a; }
        if (side[sbd].eb == eb) { side[sbd].eb = ea; side[sbd].b = a; }

        circles(ea);
        circles(eb);
    }
    /*-swap_side---------------------------------------------------------------*/

    void circles(int e)
    /*---------------------------------------------------+
    |  This function calculates radii of inscribed and   |
    |  circumscribed circle for a given element (int e)  |
    +---------------------------------------------------*/
    {
        float x, y, xi, yi, xj, yj, xk, yk, xij, yij, xjk, yjk, num, den;
        float si, sj, sk, O;

        xi = node[elem[e].i].x; yi = node[elem[e].i].y;
        xj = node[elem[e].j].x; yj = node[elem[e].j].y;
        xk = node[elem[e].k].x; yk = node[elem[e].k].y;

        xij = 0.5f * (xi + xj); yij = 0.5f * (yi + yj);
        xjk = 0.5f * (xj + xk); yjk = 0.5f * (yj + yk);

        num = (xij - xjk) * (xj - xi) + (yij - yjk) * (yj - yi);
        den = (xj - xi) * (yk - yj) - (xk - xj) * (yj - yi);

        if (den > 0)
        {
            elem[e].xv = x = xjk + num / den * (yk - yj);
            elem[e].yv = y = yjk - num / den * (xk - xj);

            elem[e].R = Mathf.Sqrt((xi - x) * (xi - x) + (yi - y) * (yi - y));
        }

        si = side[elem[e].si].s;
        sj = side[elem[e].sj].s;
        sk = side[elem[e].sk].s;
        O = si + sj + sk;
        elem[e].Det = xi * (yj - yk) - xj * (yi - yk) + xk * (yi - yj);

        elem[e].xin = (xi * si + xj * sj + xk * sk) / O;
        elem[e].yin = (yi * si + yj * sj + yk * sk) / O;

        elem[e].r = elem[e].Det / O;
    }

    void spacing(int e, int n)
    /*----------------------------------------------------------------+
    |  This function calculates the value of the spacing function in  |
    |  a new node 'n' which is inserted in element 'e' by a linear    |
    |  approximation from the values of the spacing function in the   |
    |  elements nodes.                                                |
    +----------------------------------------------------------------*/
    {
        float dxji, dxki, dyji, dyki, dx_i, dy_i, det, a, b;

        dxji = node[elem[e].j].x - node[elem[e].i].x;
        dyji = node[elem[e].j].y - node[elem[e].i].y;
        dxki = node[elem[e].k].x - node[elem[e].i].x;
        dyki = node[elem[e].k].y - node[elem[e].i].y;
        dx_i = node[n].x - node[elem[e].i].x;
        dy_i = node[n].y - node[elem[e].i].y;

        det = dxji * dyki - dxki * dyji;

        a = (+dyki * dx_i - dxki * dy_i) / det;
        b = (-dyji * dx_i + dxji * dy_i) / det;

        node[n].F = node[elem[e].i].F + 
            a * (node[elem[e].j].F - node[elem[e].i].F) +  
            b * (node[elem[e].k].F - node[elem[e].i].F);
    }

    int insert_node(float x, float y, int spac,
     int prev_n, int prev_s_mark, int mark, int next_s_mark, int next_n) 
    {
        int i, j, k, en, n, e, ei, ej, ek, s, si, sj, sk;
        float sx, sy;

        Nn++;          /* one new node */

        node[Nn - 1].x = x;
        node[Nn - 1].y = y;
        node[Nn - 1].mark = mark;

        /* find the element which contains new node */
        e = in_elem(node[Nn - 1]);

        /* calculate the spacing function in the new node */
        if (spac == ON)
            spacing(e, Nn - 1);

        i = elem[e].i; j = elem[e].j; k = elem[e].k;
        ei = elem[e].ei; ej = elem[e].ej; ek = elem[e].ek;
        si = elem[e].si; sj = elem[e].sj; sk = elem[e].sk;

        Ne += 2;
        Ns += 3;

        /*---------------+
        |  new elements  |
        +---------------*/
        elem[Ne - 2].i = Nn - 1; elem[Ne - 2].j = k; elem[Ne - 2].k = i;
        elem[Ne - 1].i = Nn - 1; elem[Ne - 1].j = i; elem[Ne - 1].k = j;

        elem[Ne - 2].ei = ej; elem[Ne - 2].ej = Ne - 1; elem[Ne - 2].ek = e;
        elem[Ne - 1].ei = ek; elem[Ne - 1].ej = e; elem[Ne - 1].ek = Ne - 2;

        elem[Ne - 2].si = sj; elem[Ne - 2].sj = Ns - 2; elem[Ne - 2].sk = Ns - 3;
        elem[Ne - 1].si = sk; elem[Ne - 1].sj = Ns - 1; elem[Ne - 1].sk = Ns - 2;

        /*------------+ 
        |  new sides  |
        +------------*/
        side[Ns - 3].c = k; side[Ns - 3].d = Nn - 1;     /* c-d */
        side[Ns - 3].a = j; side[Ns - 3].b = i;        /* a-b */
        side[Ns - 3].ea = e; side[Ns - 3].eb = Ne - 2;

        side[Ns - 2].c = i; side[Ns - 2].d = Nn - 1;     /* c-d */
        side[Ns - 2].a = k; side[Ns - 2].b = j;        /* a-b */
        side[Ns - 2].ea = Ne - 2; side[Ns - 2].eb = Ne - 1;

        side[Ns - 1].c = j; side[Ns - 1].d = Nn - 1;     /* c-d */
        side[Ns - 1].a = i; side[Ns - 1].b = k;        /* a-b */
        side[Ns - 1].ea = Ne - 1; side[Ns - 1].eb = e;

        for (s = 1; s <= 3; s++)
        {
            sx = node[side[Ns - s].c].x - node[side[Ns - s].d].x;
            sy = node[side[Ns - s].c].y - node[side[Ns - s].d].y;
            side[Ns - s].s = Mathf.Sqrt(sx * sx + sy * sy);
        }

        elem[e].i = Nn - 1;
        elem[e].ej = Ne - 2;
        elem[e].ek = Ne - 1;
        elem[e].sj = Ns - 3;
        elem[e].sk = Ns - 1;

        if (side[si].a == i) { side[si].a = Nn - 1; side[si].ea = e; }
        if (side[si].b == i) { side[si].b = Nn - 1; side[si].eb = e; }

        if (side[sj].a == j) { side[sj].a = Nn - 1; side[sj].ea = Ne - 2; }
        if (side[sj].b == j) { side[sj].b = Nn - 1; side[sj].eb = Ne - 2; }

        if (side[sk].a == k) { side[sk].a = Nn - 1; side[sk].ea = Ne - 1; }
        if (side[sk].b == k) { side[sk].b = Nn - 1; side[sk].eb = Ne - 1; }

        if (ej != -1)
        {
            if (elem[ej].ei == e) { elem[ej].ei = Ne - 2; }
            if (elem[ej].ej == e) { elem[ej].ej = Ne - 2; }
            if (elem[ej].ek == e) { elem[ej].ek = Ne - 2; }
        }

        if (ek != -1)
        {
            if (elem[ek].ei == e) { elem[ek].ei = Ne - 1; }
            if (elem[ek].ej == e) { elem[ek].ej = Ne - 1; }
            if (elem[ek].ek == e) { elem[ek].ek = Ne - 1; }
        }

        /* Find circumenters for two new elements, 
           and for the one who's segment has changed */
        circles(e);
        circles(Ne - 2);
        circles(Ne - 1);

        bowyer(Nn - 1, spac);

        /*-------------------------------------------------+
        |  NEW ! Insert boundary conditions for the sides  |
        +-------------------------------------------------*/
        for (s = 3; s < Ns; s++)
        {
            if (side[s].c == prev_n && side[s].d == Nn - 1) side[s].mark = prev_s_mark;
            if (side[s].d == prev_n && side[s].c == Nn - 1) side[s].mark = prev_s_mark;
            if (side[s].c == next_n && side[s].d == Nn - 1) side[s].mark = next_s_mark;
            if (side[s].d == next_n && side[s].c == Nn - 1) side[s].mark = next_s_mark;
        }

        return e;
    }

    void erase()
    {
        int s, n, e;

        int a, b, c, d, ea, eb;

        /*--------------------------+
        |                           |
        |  Negative area check for  |
        |  elimination of elements  |
        |                           |
        +--------------------------*/
        for (e = 0; e < Ne; e++)
            if ((node[elem[e].i].chain == node[elem[e].j].chain) &&
                (node[elem[e].j].chain == node[elem[e].k].chain) &&
                (chain[node[elem[e].i].chain].type == CLOSED))
            {
                a = Mathf.Min(Mathf.Min(elem[e].i, elem[e].j), elem[e].k);
                c = Mathf.Max(Mathf.Max(elem[e].i, elem[e].j), elem[e].k);
                b = elem[e].i + elem[e].j + elem[e].k - a - c;

                if (a < 3)
                    elem[e].mark = OFF;

                else if (area(node[a], node[b], node[c]) < 0.0f)
                    elem[e].mark = OFF;
            }

        for (e = 0; e < Ne; e++)
        {
            if (elem[elem[e].ei].mark == OFF) elem[e].ei = OFF;
            if (elem[elem[e].ej].mark == OFF) elem[e].ej = OFF;
            if (elem[elem[e].ek].mark == OFF) elem[e].ek = OFF;
        }

        /*-----------------------+
        |                        |
        |  Elimination of sides  |
        |                        |
        +-----------------------*/
        for (s = 0; s < 3; s++)
            side[s].mark = OFF;

        for (s = 3; s < Ns; s++)
            if ((elem[side[s].ea].mark == OFF) && (elem[side[s].eb].mark == OFF))
                side[s].mark = OFF;

        for (s = 3; s < Ns; s++)
            if (side[s].mark != OFF)
            {
                if (elem[side[s].ea].mark == OFF) { side[s].ea = OFF; side[s].a = OFF; }
                if (elem[side[s].eb].mark == OFF) { side[s].eb = OFF; side[s].b = OFF; }
            }

        /*-----------------------+
        |                        |
        |  Elimination of nodes  |
        |                        |
        +-----------------------*/
        for (n = 0; n < 3; n++)
            node[n].mark = OFF;

    }

    void diamond()
    {
        int ea, eb, s;
        int eac = 0;
        int ead = 0;
        int ebc = 0;
        int ebd = 0;

        for (s = 0; s < Ns; s++)
            if (side[s].mark != OFF)
            {
                ea = side[s].ea;
                eb = side[s].eb;

                if (elem[ea].ei == eb) { ead = elem[ea].ej; eac = elem[ea].ek; }
                if (elem[ea].ej == eb) { ead = elem[ea].ek; eac = elem[ea].ei; }
                if (elem[ea].ek == eb) { ead = elem[ea].ei; eac = elem[ea].ej; }
                if (elem[eb].ei == ea) { ebc = elem[eb].ej; ebd = elem[eb].ek; }
                if (elem[eb].ej == ea) { ebc = elem[eb].ek; ebd = elem[eb].ei; }
                if (elem[eb].ek == ea) { ebc = elem[eb].ei; ebd = elem[eb].ej; }

                if ((eac == OFF || elem[eac].state == D) &&
                (ebc == OFF || elem[ebc].state == D) &&
                (ead == OFF || elem[ead].state == D) &&
                (ebd == OFF || elem[ebd].state == D))
                {
                    elem[ea].state = D;
                    elem[eb].state = D;
                }
            }
    }

    void classify()
    /*----------------------------------------------------------+
    |  This function searches through all elements every time.  |
    |  Some optimisation will definitely bee needed             |
    |                                                           |
    |  But it also must me noted, that this function defines    |
    |  the strategy for insertion of new nodes                  |
    |                                                           |
    |  It's MUCH MUCH better when the ugliest element is found  |
    |  as one with highest ratio of R/r !!! (before it was      |
    |  element with greater R)                                  |
    +----------------------------------------------------------*/
    {
        int e, ei, ej, ek, si, sj, sk;
        float ratio = -GREAT, F;

        ugly = OFF;

        for (e = 0; e < Ne; e++)
            if (elem[e].mark != OFF)
            {
                ei = elem[e].ei; ej = elem[e].ej; ek = elem[e].ek;

                F = (node[elem[e].i].F + node[elem[e].j].F + node[elem[e].k].F) / 3.0f;

                elem[e].state = W;

                /*--------------------------+
                |  0.577 is ideal triangle  |
                +--------------------------*/
                if (elem[e].R < 0.700 * F) elem[e].state = D; /* 0.0866; 0.07 */

                /*------------------------+
                |  even this is possible  |
                +------------------------*/
                if (ei != OFF && ej != OFF && ek != OFF)
                    if (elem[ei].state == D && elem[ej].state == D && elem[ek].state == D)
                        elem[e].state = D;
            }

        /*--------------------------------------+
        |  Diamond check. Is it so important ?  |
        +--------------------------------------*/
        diamond();

        /*------------------------------------------------+
        |  First part of the trick:                       |
        |    search through the elements on the boundary  |
        +------------------------------------------------*/
        for (e = 0; e < Ne; e++)
            if (elem[e].mark != OFF && elem[e].state != D)
            {
                si = elem[e].si; sj = elem[e].sj; sk = elem[e].sk;

                if (side[si].mark != 0) elem[e].state = A;
                if (side[sj].mark != 0) elem[e].state = A;
                if (side[sk].mark != 0) elem[e].state = A;

                if (elem[e].state == A && elem[e].R / elem[e].r > ratio)
                {
                    ratio = Mathf.Max(ratio, elem[e].R / elem[e].r);
                    ugly = e;
                }
            }

        /*-------------------------------------------------+
        |  Second part of the trick:                       |
        |    if non-acceptable element on the boundary is  |
        |    found, ignore the elements inside the domain  |
        +-------------------------------------------------*/
        if (ugly == OFF)
            for (e = 0; e < Ne; e++)
                if (elem[e].mark != OFF)
                {
                    if (elem[e].state != D)
                    {
                        ei = elem[e].ei; ej = elem[e].ej; ek = elem[e].ek;

                        if (ei != OFF)
                            if (elem[ei].state == D) elem[e].state = A;

                        if (ej != OFF)
                            if (elem[ej].state == D) elem[e].state = A;

                        if (ek != OFF)
                            if (elem[ek].state == D) elem[e].state = A;

                        if (elem[e].state == A && elem[e].R / elem[e].r > ratio)
                        {
                            ratio = Mathf.Max(ratio, elem[e].R / elem[e].r);
                            ugly = e;
                        }
                    }
                }

    }

    void new_node()
    /*---------------------------------------------------+
    |  This function is very important.                  |
    |  It determines the position of the inserted node.  |
    +---------------------------------------------------*/
    {
        int s = OFF,  e;
        int n = 0;
        float xM, yM, xCa, yCa, p, px, py, q, qx, qy, rhoM, rho_M, d;

        nod Ca = new nod();

    /*-------------------------------------------------------------------------+
    |  It's obvious that elements which are near the boundary, will come into  |
    |  play first.                                                             |
    |                                                                          |
    |  However, some attention has to be payed for the case when two accepted  |
    |  elements surround the ugly one                                          |
    |                                                                          |
    |  What if new points falls outside the domain                             |
    +-------------------------------------------------------------------------*/
     if(elem[elem[ugly].ei].state==D)    {s=elem[ugly].si; n=elem[ugly].i;}
     if(elem[elem[ugly].ej].state==D)    {s=elem[ugly].sj; n=elem[ugly].j;}
     if(elem[elem[ugly].ek].state==D)    {s=elem[ugly].sk; n=elem[ugly].k;}
     if(side[elem[ugly].si].mark > 0)    {s=elem[ugly].si; n=elem[ugly].i;}
     if(side[elem[ugly].sj].mark > 0)    {s=elem[ugly].sj; n=elem[ugly].j;}
     if(side[elem[ugly].sk].mark > 0)    {s=elem[ugly].sk; n=elem[ugly].k;}
     if(s==OFF) return;

     xM  = 0.5f*(node[side[s].c].x + node[side[s].d].x);
     yM  = 0.5f*(node[side[s].c].y + node[side[s].d].y);

     Ca.x = elem[ugly].xv;
     Ca.y = elem[ugly].yv;

     p  = 0.5f* side[s].s;    /* not checked */

     qx = Ca.x-xM;
     qy = Ca.y-yM;
     q  = Mathf.Sqrt(qx* qx+qy* qy);

     rhoM = 0.577f *  0.5f*(node[side[s].c].F + node[side[s].d].F);

     rho_M = Mathf.Min(Mathf.Max( rhoM, p), 0.5f*(p* p+q* q)/q );

     if(rho_M<p)
            d =rho_M;
     else
            d =rho_M+Mathf.Sqrt(rho_M* rho_M-p* p); 

    /*---------------------------------------------------------------------+
    |  The following line check can the new point fall outside the domain. |
    |  However, I can't remember how it works, but I believe that it is    |
    |  still a weak point of the code, particulary when there are lines    |
    |  inside the domain                                                   |
    +---------------------------------------------------------------------*/

     if( area(node[side[s].c], node[side[s].d], Ca) *
         area(node[side[s].c], node[side[s].d], node[n]) > 0.0f )
       insert_node(xM + d* qx/q, yM + d* qy/q, ON, OFF, 0, 0, 0, OFF);
    /*
     else
      {
       node[n].x = xM - d*qx/q;
       node[n].y = yM - d*qy/q;
       node[n].mark=6;   
       for(e=0; e<Ne; e++) 
         if(elem[e].i==n || elem[e].j==n || elem[e].k==n)
           circles(e);
      }
    */
     return;
    }

    void neighbours()
    /*--------------------------------------------------------------+
    |  Counting the elements which surround each node.              |
    |  It is important for the two functions: 'relax' and 'smooth'  |
    +--------------------------------------------------------------*/
    {
        int s;

        for (s = 0; s < Ns; s++)
            if (side[s].mark == 0)
            {
                if (node[side[s].c].mark == 0)
                    node[side[s].c].Nne++;

                if (node[side[s].d].mark == 0)
                    node[side[s].d].Nne++;
            }
    }

    void materials()
    {
        int e, c, over; //iter, , s
        int mater = 0;
        int ei, ej, ek, si, sj, sk;

        for (e = 0; e < Ne; e++)
            if (elem[e].mark != OFF)
                elem[e].material = OFF;

        for (c = 0; c < Nc; c++)
        {
            if (point[c].inserted == 0)
            {
                elem[in_elem(point[c])].material = point[c].mark;
                mater = ON;
            }
        }

        if (mater == ON)
        {
            for (;;)
            {
                over = ON;

                for (e = 0; e < Ne; e++)
                    if (elem[e].mark != OFF && elem[e].material == OFF)
                    {
                        ei = elem[e].ei;
                        ej = elem[e].ej;
                        ek = elem[e].ek;

                        si = elem[e].si;
                        sj = elem[e].sj;
                        sk = elem[e].sk;


                        if (ei != OFF)
                            if (elem[ei].material != OFF && side[si].mark == 0)
                            {
                                elem[e].material = elem[ei].material;
                                over = OFF;
                            }

                        if (ej != OFF)
                            if (elem[ej].material != OFF && side[sj].mark == 0)
                            {
                                elem[e].material = elem[ej].material;
                                over = OFF;
                            }

                        if (ek != OFF)
                            if (elem[ek].material != OFF && side[sk].mark == 0)
                            {
                                elem[e].material = elem[ek].material;
                                over = OFF;
                            }

                    }

                if (over == ON) break;

            } /* for(iter) */

        }
    }

    void relax()
    {
        int s, T, E;

        for (T = 6; T >= 3; T--)
        { 
            for (s = 0; s < Ns; s++)
            {
                if (side[s].mark == 0)
                {

                    if ((node[side[s].a].mark == 0) &&
                    (node[side[s].b].mark == 0) &&
                    (node[side[s].c].mark == 0) &&
                    (node[side[s].d].mark == 0))
                    {
                        E = node[side[s].c].Nne + node[side[s].d].Nne
                        - node[side[s].a].Nne - node[side[s].b].Nne;

                        if (E == T)
                        {
                            node[side[s].a].Nne++; node[side[s].b].Nne++;
                            node[side[s].c].Nne--; node[side[s].d].Nne--;
                            swap_side(s);
                        }
                    }
                }
	        }
        }
    }

    int smooth()
    {
        int it, s, n, e;

        for (it = 0; it < 10; it++)
        {
            for (s = 0; s < Ns; s++)
            {
	            if (side[s].mark == 0)
	            {
	                if (node[side[s].c].mark == 0)
	                {
	                    node[side[s].c].sumx += node[side[s].d].x;
	                    node[side[s].c].sumy += node[side[s].d].y;
	                }
	
	                if (node[side[s].d].mark == 0)
	                {
	                    node[side[s].d].sumx += node[side[s].c].x;
	                    node[side[s].d].sumy += node[side[s].c].y;
	                }
	            }
            }

            for (n = 0; n < Nn; n++) { 
                if (node[n].mark == 0)
                {
                    node[n].x = node[n].sumx / node[n].Nne; node[n].sumx = 0.0f;
                    node[n].y = node[n].sumy / node[n].Nne; node[n].sumy = 0.0f;
                }
            }
        }

        for (e = 0; e < Ne; e++)
            if (elem[e].mark != OFF)
                circles(e);

        return 0;
    }

    void renum()
    {
        int n, o, s, e, e2, c, d, i, j, k;
        int new_elem = 0, new_node = 0, new_side = 0, next_e, next_s, lowest;

        for (n = 0; n < Nn; n++) node[n].new_numb = OFF;
        for (e = 0; e < Ne; e++) elem[e].new_numb = OFF;
        for (s = 0; s < Ns; s++) side[s].new_numb = OFF;

        /*-------------------------------+
        |  Searching the first element.  |
        |  It is the first which is ON   |
        +-------------------------------*/
        for (e = 0; e < Ne; e++)
            if (elem[e].mark != OFF)
                break;

        /*----------------------------------------------------------+
        |  Assigning numbers 0 and 1 to the nodes of first element  |
        +----------------------------------------------------------*/
        node[elem[e].i].new_numb = new_node; new_node++;
        node[elem[e].j].new_numb = new_node; new_node++;

        /*%%%%%%%%%%%%%%%%%%%%%%%%%
        %                         %
        %  Renumeration of nodes  %
        %                         % 
        %%%%%%%%%%%%%%%%%%%%%%%%%*/
        do
        {
            lowest = Nn + Nn;
            next_e = OFF;

            for (e = 0; e < Ne; e++)
                if (elem[e].mark != OFF && elem[e].new_numb == OFF)
                {
                    i = node[elem[e].i].new_numb;
                    j = node[elem[e].j].new_numb;
                    k = node[elem[e].k].new_numb;

                    if (i + j + k + 2 == Mathf.Abs(i) + Mathf.Abs(j) + Mathf.Abs(k))
                    {
                        if ((i == OFF) && (j + k) < lowest) { next_e = e; lowest = j + k; }
                        if ((j == OFF) && (i + k) < lowest) { next_e = e; lowest = i + k; }
                        if ((k == OFF) && (i + j) < lowest) { next_e = e; lowest = i + j; }
                    }
                }

            if (next_e != OFF)
            {
                i = node[elem[next_e].i].new_numb;
                j = node[elem[next_e].j].new_numb;
                k = node[elem[next_e].k].new_numb;

                /*----------------------------------+
                |  Assign a new number to the node  |
                +----------------------------------*/
                if (i == OFF) { node[elem[next_e].i].new_numb = new_node; new_node++; }
                if (j == OFF) { node[elem[next_e].j].new_numb = new_node; new_node++; }
                if (k == OFF) { node[elem[next_e].k].new_numb = new_node; new_node++; }
            }
        }
        while (next_e != OFF);

        /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                             %
        %  Renumeration of triangles  %
        %                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
        do
        {
            lowest = Nn + Nn + Nn;
            next_e = OFF;

            for (e = 0; e < Ne; e++)
                if (elem[e].mark != OFF && elem[e].new_numb == OFF)
                {
                    i = node[elem[e].i].new_numb;
                    j = node[elem[e].j].new_numb;
                    k = node[elem[e].k].new_numb;

                    if ((i + j + k) < lowest)
                    {
                        next_e = e;
                        lowest = i + j + k;
                    }
                }

            if (next_e != OFF)
            {
                elem[next_e].new_numb = new_elem; new_elem++;
            }
        }
        while (next_e != OFF);



        /*%%%%%%%%%%%%%%%%%%%%%%%%%
        %                         %
        %  Renumeration of sides  %
        %                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%*/
        do
        {
            lowest = Nn + Nn;
            next_s = OFF;

            for (s = 0; s < Ns; s++)
                if (side[s].mark != OFF && side[s].new_numb == OFF)
                {
                    c = node[side[s].c].new_numb;
                    d = node[side[s].d].new_numb;

                    if ((c + d) < lowest)
                    {
                        lowest = c + d; next_s = s;
                    }
                }

            if (next_s != OFF)
            {
                side[next_s].new_numb = new_side;
                new_side++;
            }

        }
        while (next_s != OFF);

    }


}




