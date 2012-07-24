#ifndef __PYTOAST_REGMGR_H
#define __PYTOAST_REGMGR_H

/**
 * \brief Manages a list of persistent regularisation objects
 */
class RegularisationManager {
public:
    RegularisationManager ();
    int Add (Regularisation *reg);
    Regularisation *Get (int idx) const;
    bool Delete (int idx);
    void Clear ();

 private:
    Regularisation **list;
    int nlist;
    int nreg;
};

#endif // !__PYTOAST_REGMGR_H
