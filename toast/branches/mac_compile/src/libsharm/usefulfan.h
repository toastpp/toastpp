#ifndef __USEFULFAN_H
#define __USEFULFAN_H


//***************************************************************************
//***************************************************************************
//*   Functions Useful for the optimisation                                 *
//***************************************************************************
//***************************************************************************
//***************************************************************************



double atan3(double hmit, double syn);// The ArcTan

RVector cross_product(RVector edge1,RVector edge2);//Cross Product for RVectors

double Sij(double x1, double y1,double z1,double x2,double y2,double z2);// the Side of the spherical quadrilaterals

double  SijV(RVector  x1,RVector x2);// the Side of the spherical quadrilaterals
double NormS(double a, double b, double c);//norm using three doubles
double NormV(RVector x);// The Magnitude of the RVector

double  Sin_Side(RVector edge1,RVector edge2);// The sine of the Side 

double SinLDAB(RVector VecD, RVector VecA, RVector VecB);//The Sine of the angle DAB

double CosLDAB(RVector VecD, RVector VecA, RVector VecB);//The cossine of the angle DAB

double  LDAB (RVector VecD, RVector VecA, RVector VecB);//The angle DAB

double Area( RVector VecA, RVector VecB, RVector  VecC, RVector VecD);//The area of the spherical quadrilateral

double DAreaDX( RVector a, RVector b, RVector  c, RVector d);// the partial derivative of the area function in respect to DXa

double DAreaDY( RVector a, RVector b, RVector  c, RVector d);// the partial derivative of the area function in respect to DYa

double DAreaDZ( RVector a, RVector b, RVector  c, RVector d);// the partial derivative of the area function in respect to DZa


double DAx(RVector VecA, RVector VecB, RVector VecC);// The partial derivatives of the angle ABC functions
double DAy(RVector VecA, RVector VecB, RVector VecC);
double DAz(RVector VecA, RVector VecB, RVector VecC);
double DBx(RVector VecA, RVector VecB, RVector VecC);
double DBy(RVector VecA, RVector VecB, RVector VecC);
double DBz(RVector VecA, RVector VecB, RVector VecC);
double DCx(RVector VecA, RVector VecB, RVector VecC);
double DCy(RVector VecA, RVector VecB, RVector VecC);
double DCz(RVector VecA, RVector VecB, RVector VecC);

double dabs(double in);//abs for doubles


#endif
