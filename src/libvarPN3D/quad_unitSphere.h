/* Weights for the quadrature scheme (which is accurate till 17th order) required for boundary integrals
*/
double wts17[] = {.0038282704949371616, .0038282704949371616, .0038282704949371616, .0038282704949371616, .0038282704949371616, .0038282704949371616, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283}; 

/* Points for the quadrature scheme (which is accurate till 17th order) required for boundary integrals
*/
double pts17[110][3] = {{0, 0, 1.0}, {0, 1.0, 0}, {1.0, 0, 0}, {0, 0, -1.0}, {0, -1.0, 0}, {-1.0, 0, 0}, 

{.57735026918962576, .57735026918962576, .57735026918962576}, {.57735026918962576, .57735026918962576, -0.57735026918962576}, {.57735026918962576, -0.57735026918962576, .57735026918962576}, {-0.57735026918962576, .57735026918962576, .57735026918962576}, {-0.57735026918962576, -0.57735026918962576, .57735026918962576}, {-0.57735026918962576, .57735026918962576, -0.57735026918962576}, {.57735026918962576, -0.57735026918962576, -0.57735026918962576}, {-0.57735026918962576, -0.57735026918962576, -0.57735026918962576}, 

{0.18511563534473617, 0.18511563534473617, 0.9651240350865941}, {0.18511563534473617, 0.18511563534473617, -0.9651240350865941}, {0.18511563534473617, -0.18511563534473617, 0.9651240350865941}, {-0.18511563534473617, 0.18511563534473617, 0.9651240350865941}, {-0.18511563534473617, -0.18511563534473617, 0.9651240350865941}, {-0.18511563534473617, 0.18511563534473617, -0.9651240350865941}, {0.18511563534473617, -0.18511563534473617, -0.9651240350865941}, {-0.18511563534473617, -0.18511563534473617, -0.9651240350865941}, 

{0.18511563534473617, 0.9651240350865941, 0.18511563534473617}, {0.18511563534473617, 0.9651240350865941, -0.18511563534473617}, {0.18511563534473617, -0.9651240350865941, 0.18511563534473617}, {-0.18511563534473617, 0.9651240350865941, 0.18511563534473617}, {-0.18511563534473617, -0.9651240350865941, 0.18511563534473617}, {-0.18511563534473617, 0.9651240350865941, -0.18511563534473617}, {0.18511563534473617, -0.9651240350865941, -0.18511563534473617}, {-0.18511563534473617, -0.9651240350865941, -0.18511563534473617}, 

{0.9651240350865941, 0.18511563534473617, 0.18511563534473617}, {0.9651240350865941, 0.18511563534473617, -0.18511563534473617}, {0.9651240350865941, -0.18511563534473617, 0.18511563534473617}, {-0.9651240350865941, 0.18511563534473617, 0.18511563534473617}, {-0.9651240350865941, -0.18511563534473617, 0.18511563534473617}, {-0.9651240350865941, 0.18511563534473617, -0.18511563534473617}, {0.9651240350865941, -0.18511563534473617, -0.18511563534473617}, {-0.9651240350865941, -0.18511563534473617, -0.18511563534473617}, 

{0.39568947305594191, 0.39568947305594191, 0.8287699812525922}, {0.39568947305594191, 0.39568947305594191, -0.8287699812525922}, {0.39568947305594191, -0.39568947305594191, 0.8287699812525922}, {-0.39568947305594191, 0.39568947305594191, 0.8287699812525922}, {-0.39568947305594191, -0.39568947305594191, 0.8287699812525922}, {-0.39568947305594191, 0.39568947305594191, -0.8287699812525922}, {0.39568947305594191, -0.39568947305594191, -0.8287699812525922}, {-0.39568947305594191, -0.39568947305594191, -0.8287699812525922}, 

{0.39568947305594191, 0.8287699812525922, 0.39568947305594191}, {0.39568947305594191, 0.8287699812525922, -0.39568947305594191}, {0.39568947305594191, -0.8287699812525922, 0.39568947305594191}, {-0.39568947305594191, 0.8287699812525922, 0.39568947305594191}, {-0.39568947305594191, -0.8287699812525922, 0.39568947305594191}, {-0.39568947305594191, 0.8287699812525922, -0.39568947305594191}, {0.39568947305594191, -0.8287699812525922, -0.39568947305594191}, {-0.39568947305594191, -0.8287699812525922, -0.39568947305594191}, 

{0.8287699812525922, 0.39568947305594191, 0.39568947305594191}, {0.8287699812525922, 0.39568947305594191, -0.39568947305594191}, {0.8287699812525922, -0.39568947305594191, 0.39568947305594191}, {-0.8287699812525922, 0.39568947305594191, 0.39568947305594191}, {-0.8287699812525922, -0.39568947305594191, 0.39568947305594191}, {-0.8287699812525922, 0.39568947305594191, -0.39568947305594191}, {0.8287699812525922, -0.39568947305594191, -0.39568947305594191}, {-0.8287699812525922, -0.39568947305594191, -0.39568947305594191}, 

{0.69042104838229218, 0.69042104838229218, 0.21595729184584883}, {0.69042104838229218, 0.69042104838229218, -0.21595729184584883}, {0.69042104838229218, -0.69042104838229218, 0.21595729184584883}, {-0.69042104838229218, 0.69042104838229218, 0.21595729184584883}, {-0.69042104838229218, -0.69042104838229218, 0.21595729184584883}, {-0.69042104838229218, 0.69042104838229218, -0.21595729184584883}, {0.69042104838229218, -0.69042104838229218, -0.21595729184584883}, {-0.69042104838229218, -0.69042104838229218, -0.21595729184584883},

{0.69042104838229218, 0.21595729184584883, 0.69042104838229218}, {0.69042104838229218, 0.21595729184584883, -0.69042104838229218}, {0.69042104838229218, -0.21595729184584883, 0.69042104838229218}, {-0.69042104838229218, 0.21595729184584883, 0.69042104838229218}, {-0.69042104838229218, -0.21595729184584883, 0.69042104838229218}, {-0.69042104838229218, 0.21595729184584883, -0.69042104838229218}, {0.69042104838229218, -0.21595729184584883, -0.69042104838229218}, {-0.69042104838229218, -0.21595729184584883, -0.69042104838229218},

{0.21595729184584883, 0.69042104838229218, 0.69042104838229218}, {0.21595729184584883, 0.69042104838229218, -0.69042104838229218}, {0.21595729184584883, -0.69042104838229218, 0.69042104838229218}, {-0.21595729184584883, 0.69042104838229218, 0.69042104838229218}, {-0.21595729184584883, -0.69042104838229218, 0.69042104838229218}, {-0.21595729184584883, 0.69042104838229218, -0.69042104838229218}, {0.21595729184584883, -0.69042104838229218, -0.69042104838229218}, {-0.21595729184584883, -0.69042104838229218, -0.69042104838229218},

{0.47836902881215020, 0, 0.8781589106040662}, {0.47836902881215020, 0, -0.8781589106040662}, {-0.47836902881215020, 0, 0.8781589106040662}, {-0.47836902881215020, 0, -0.8781589106040662}, {0.8781589106040662, 0, 0.47836902881215020}, {0.8781589106040662, 0, -0.47836902881215020}, {-0.8781589106040662, 0, 0.47836902881215020}, {-0.8781589106040662, 0, -0.47836902881215020}, 

{0, 0.8781589106040662, 0.47836902881215020}, {0, 0.8781589106040662, -0.47836902881215020}, {0, -0.8781589106040662, 0.47836902881215020}, {0, -0.8781589106040662, -0.47836902881215020}, {0, 0.47836902881215020, 0.8781589106040662}, {0, 0.47836902881215020, -0.8781589106040662}, {0, -0.47836902881215020, 0.8781589106040662}, {0, -0.47836902881215020, -0.8781589106040662}, 

{0.47836902881215020, 0.8781589106040662, 0}, {0.47836902881215020, -0.8781589106040662, 0}, {-0.47836902881215020, 0.8781589106040662, 0}, {-0.47836902881215020, -0.8781589106040662, 0}, {0.8781589106040662, 0.47836902881215020, 0}, {0.8781589106040662, -0.47836902881215020, 0}, {-0.8781589106040662, 0.47836902881215020, 0}, {-0.8781589106040662, -0.47836902881215020, 0}};
