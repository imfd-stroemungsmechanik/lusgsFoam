// Read settings for LUSGS

// Maximum number of iterations
label maxIter(readLabel(mesh.solutionDict().subDict("LUSGS").lookup("maxIter")));
// Relative solver tolerance
scalar relTol(readScalar(mesh.solutionDict().subDict("LUSGS").lookup("relTol")));

// Initial and current residuals
scalar initRho(0.0);
scalar initRhoU(0.0);
scalar initRhoE(0.0);
scalar curRho(0.0);
scalar curRhoU(0.0);
scalar curRhoE(0.0);
