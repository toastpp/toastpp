#ifndef __PYTOAST_RASTERMGR_H
#define __PYTOAST_RASTERMGR_H

/**
 * \brief Manages a list of persistent active basis mappers
 */
class RasterManager {
public:
    RasterManager ();
    int Add (Raster *raster);
    Raster *Get (int idx) const;
    bool Delete (int idx);
    void Clear ();

private:
    Raster **list;
    int nlist;
    int nraster;
};

#endif // !__PYTOAST_RASTERMGR_H
