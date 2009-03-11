algebraic3d

solid slab = plane (-50, -50, 0; 0, 0, -1)
	and plane (-50, -50, 0; 0, -1, 0)
	and plane (-50, -50, 0; -1, 0, 0)
	and plane (50, 50, 20; 0, 0, 1)
	and plane (50, 50, 20; 0, 1, 0)
	and plane (50, 50, 20; 1, 0, 0);

solid cyl = cylinder (-50, 0, 10; 50, 0, 10; 7 )
	and plane (-50, 0, 10; -1, 0, 0)
	and plane ( 50, 0, 10; 1, 0, 0);

solid slabandhole = slab and not cyl;

tlo slabandhole -transparent -col=[0,0,1]; 
tlo cyl -col=[1,0,0];
