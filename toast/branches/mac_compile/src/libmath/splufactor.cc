template<class MT>
int spOrderAndFactor (TCompRowMatrix<MT> &Matrix, TVector<MT> &rhs,
    double RelThreshold, double AbsThreshold, bool DiagPivoting)
{
    int Step, Size;
    MT *pPivot;
    Size = Matrix.nRows();

    if (!Matrix.NeedsOrdering) {
        for (Step = 0; Step < Size; Step++) {
}
