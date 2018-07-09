#############################################################################
##  orders package - "orders.gd"
##  Copyright 2018 Florian Eisele
#############################################################################

## Declarations:

ZpOrderFamily := NewFamily("ZpOrderFamily");
DeclareCategory("IsZpOrder", IsObject); # Z_(p)-orders in semisimple (!) Q_(p)-algebras

#!! IsZpOrderMatrixRep is deprecated. Use IsZpOrderMultiMatrixRep instead!
DeclareRepresentation("IsZpOrderMatrixRep", IsComponentObjectRep, [ "p", "n", "gens" ]);
	# In addition to that, the following components of IsZpOrderMatrixRep may be defined:
	# "simple":                List of the simple Lambda-modules
	# "simple_multiplicities": List of their multiplicities (should always be defined when "simple" is defined)
     # "Qp_end_basis":          Qp-basis of the K(x)Lambda-endomorphism-ring of the representation given by "gens"
     # Qp_end_basis_smallmats:  A lower dimensional representation of Qp_end_basis

DeclareRepresentation("IsZpOrderMultiMatrixRep", IsComponentObjectRep, [ "p", "nvec", "gens" ]);
	# gens is supposed to be a list of matrices, duch that List(gens, x->x[i]) are the representations
	# on the irreducible lattices
	# In addition to that, the following can be defined:
	# "simple":                 List of the simple Lambda-modules
	# "simple_multiplicities":  The decomposition matrix
	# "Qp_end_bases":           As is IsZpOrderMatrixRep, but here it is a list of bases, one for each Wedderburn-
	#                           component
	# "Qp_end_bases_smallmats": Lower-dimensional representations of the above.
	# "component_names":        List of names of the Wedderburn-components
	# "simple_names":           List of names of the simple modules
	# ZpSn and the like also assigns values to:
	# "eval_map":               List of functions from a group G into Gl(n,Z), corresponding to the
	#                           representations on IrreducibleLattices(...). This is good for faster condensation.
	# Functions to determine projectives also assign:
	# "pim":                    The porjective covers of the modules in "simple".
	# BasicOrder also assigns values to the following:
	# "idempot":                A full set of primitive idempotents (in the same arder as "simple"!)
	# "blockmap":               List with entries of the form [[i,j], pos], meaning that the positions (in gens)
	#                           of the basis elements belonging to Hom(P(S_i), P(S_j))^tr can be found in blockidx[pos].
	# "blockidx":               List with entries of the form [p, len]. These mean that gens{[p..p+len-1]} is a basis of
	#                           some Hom(P(S_i), P(S_j)). i and j can be retrieved from blockmap
	# "gens_are_basis":         Should be set to true.

DeclareGlobalFunction("ZpOrderByMatrices");
DeclareGlobalFunction("ZpOrderByMultiMatrices");
DeclareOperation("PrintAsFunction", [IsObject]);
DeclareOperation("SaveAsFunction", [IsString, IsObject]);
DeclareOperation("SaveAsRecord", [IsString, IsZpOrder]);
DeclareOperation("ReadOrder", [IsString]);
DeclareOperation("ReadModule", [IsString, IsZpOrder]);

DeclareAttribute("DecompositionMatrix", IsZpOrder);
DeclareOperation("SimpleModules", [IsZpOrder]);
DeclareOperation("IrreducibleLattices", [IsZpOrder]);
DeclareOperation("ProjectiveIndecomposableLattices", [IsZpOrder]);

DeclareOperation("NameSimpleModulesByDims", [IsZpOrder]);
DeclareOperation("NameWedderburnComponentsByDims", [IsZpOrder]);

DeclareOperation("ExtractWedderburnComponents", [IsZpOrder, IsList]);
DeclareOperation("DirectSumOfOrders", [IsList]);
DeclareOperation("BlocksOfZpOrder", [IsZpOrder]);
DeclareOperation("DirectSumOfOrders", [IsZpOrder, IsZpOrder]);

DeclareGlobalFunction("UnFlattenMultiMatrixNC");
DeclareOperation("DirectSumOfModules", [IsList]);
DeclareOperation("DirectSumOfModules", [IsObject, IsObject]);

RModuleOverZpOrderFamily := NewFamily("RModuleOverZpOrderFamily");
DeclareCategory("IsRModuleOverZpOrder", IsObject);
DeclareCategory("IsRModuleOverZpOrderModp", IsObject);
DeclareCategory("IsRLatticeOverZpOrder", IsObject);
InstallTrueMethod(IsRModuleOverZpOrder, IsRModuleOverZpOrderModp);
InstallTrueMethod(IsRModuleOverZpOrder, IsRLatticeOverZpOrder);

DeclareRepresentation("IsRModuleByRepresentationRep", IsComponentObjectRep,
                      [ "p", "n", "gens", "order" ]);
# In addition to that, the following can be defined in IsRModuleByRepresentationRep:
# embedding_into_irr_lat: List [B, [n_1, n_2, n_3, ...]] where B is an embedding in
#                         DirectSumOfModules(List([n_1, ..], k -> IrreducibleModules[k]))
# Qp_end_basis:           Basis of the Qp(x)Lambda-endomorphism-ring of the module (only if it is a lattice)
# Qp_end_basis_smallmats: Lower-dimensional representation of Qp_end_basis

DeclareRepresentation("IsZeroRModuleRep", IsComponentObjectRep, [ "p", "order" ]); # This represents a zero module.
InstallTrueMethod(IsRModuleOverZpOrderModp, IsZeroRModuleRep);
InstallTrueMethod(IsRLatticeOverZpOrder, IsZeroRModuleRep);

DeclareGlobalFunction("ZeroRModule");
DeclareAttribute("Dimension", IsRModuleOverZpOrder);

DeclareOperation("Generators", [IsZpOrder]);
DeclareOperation("Rep", [IsRModuleOverZpOrder]);
DeclareOperation("SimpleNames", [IsZpOrder]);
DeclareOperation("ComponentNames", [IsZpOrder]);
DeclareOperation("NameSimpleModules", [IsZpOrder, IsList]);
DeclareOperation("NameWedderburnComponents", [IsZpOrder, IsList]);
DeclareOperation("InstallTrivialEndomorphismRings", [IsZpOrder]);
DeclareOperation("InstallEndomorphismRingsNC", [IsZpOrder, IsList]);
DeclareOperation("CalculateEndomorphismRingsByReynoldsNC", [IsZpOrder]);


DeclareGlobalVariable("_MAGMA_EXECUTABLE@");
DeclareGlobalFunction("SetMAGMAExecutable");
DeclareOperation("CalculateEndomorphismRingsWithMAGMA", [IsZpOrder]);


DeclareGlobalFunction("ZpSn");
DeclareGlobalFunction("ZpSnWedderburnComponentsNC");
DeclareGlobalFunction("PartitionAsString");

DeclareGlobalFunction("RModuleOverZpOrder");
DeclareOperation("ReduceModP", [IsRModuleOverZpOrder]);

DeclareOperation("MaximalSubmoduleBases", [IsRModuleOverZpOrder, IsRModuleOverZpOrderModp]);
DeclareOperation("MaximalSubmoduleBases", [IsRModuleOverZpOrder]);
DeclareOperation("MaximalSubmoduleBasesMTX", [IsRModuleOverZpOrder]);
DeclareOperation("MaximalSubmodules", [IsRModuleOverZpOrder, IsRModuleOverZpOrderModp]);
DeclareOperation("MaximalSubmodules", [IsRModuleOverZpOrder]);
DeclareOperation("Hom", [IsRModuleOverZpOrder, IsRModuleOverZpOrder]);
DeclareOperation("HomToSimpleNC", [IsRModuleOverZpOrder, IsRModuleOverZpOrderModp]);
DeclareOperation("SubmoduleByBasisNC", [IsRModuleOverZpOrder, IsObject]);
DeclareOperation("BasicOrder", [IsObject]);
DeclareOperation("RadicalOfModule", [IsRModuleOverZpOrder]);
DeclareOperation("RadicalSeries", [IsRModuleOverZpOrder, IsPosInt]);

DeclareOperation("TopEpimorphism", [IsRModuleOverZpOrder]);
DeclareOperation("LiftHomomorphismNC", [IsRModuleOverZpOrder, IsMatrix, IsRModuleOverZpOrder,
	IsMatrix, IsRModuleOverZpOrder]);
DeclareOperation("LiftHomomorphismToDirectSumNC", [IsList, IsList, IsMatrix,
  IsRModuleOverZpOrder, IsMatrix, IsRModuleOverZpOrder]);
DeclareOperation("ProjectiveCover", [IsRModuleOverZpOrder]);
DeclareOperation("HellerTranslate", [IsRModuleOverZpOrder]);
DeclareOperation("HellerTranslateModular", [IsRModuleOverZpOrder]);

DeclareOperation("ProjectiveIndecomposableLatticesForBasicOrder", [IsZpOrder]);

DeclareOperation("ProjectiveIndecomposableForBasicOrder", [IsZpOrder and IsZpOrderMultiMatrixRep, IsPosInt]);
DeclareOperation("GeneratorsForBasicOrder", [IsZpOrder and IsZpOrderMultiMatrixRep]);

DeclareGlobalFunction("LatticeAlgorithm");
DeclareGlobalFunction("GlueUpNC");
DeclareGlobalFunction("LatticesWithSimpleRadQuo");
DeclareGlobalFunction("HomForLattices");

DeclareGlobalFunction("HomByReynoldsNC");
DeclareGlobalFunction("AllLatticesCond");
DeclareGlobalFunction("AllLattices");
DeclareGlobalFunction("AllLatticesInd");
DeclareOperation("IsomorphismRModules", [IsRModuleOverZpOrder, IsRModuleOverZpOrder]);
DeclareGlobalFunction("SpinningAlgorithmNC");

DeclareOperation("CondensationData", [IsGroup, IsGroup, IsCharacter]);
DeclareOperation("CondensationData", [IsGroup, IsGroup]);
DeclareOperation("CondenseGroupRingNC", [IsZpOrder, IsList, IsRecord]);
DeclareOperation("CondenseGroupRingNC", [IsZpOrder, IsRecord]);
DeclareOperation("CondensationProperties", [IsZpOrder, IsList, IsRecord]);
DeclareOperation("CondensationProperties", [IsZpOrder, IsRecord]);
DeclareGlobalFunction("CondenseMatricesWithEvalMapNC");
DeclareGlobalFunction("CondenseMatricesNC");
DeclareGlobalFunction("CondenseTorsionRepNC");

DeclareOperation("DiagonalJoin", [IsList]);
DeclareOperation("DiagonalJoin", [IsMatrix, IsMatrix]);

DeclareOperation("GramMatrixOfTrace", [IsZpOrder, IsList]);
DeclareOperation("GramMatrixOfTrace", [IsZpOrder, IsRat]);
DeclareOperation("Valuation", [IsObject, IsInt]);
DeclareOperation("RightInverse", [IsObject]);
DeclareOperation("SelfdualSuborders", [IsZpOrder, IsList]);
DeclareOperation("AreConjugate", [IsZpOrder, IsZpOrder]);

DeclareGlobalFunction("NaturalSpechtRepresentation");
DeclareGlobalFunction("EnumerateStandardTableaux");
DeclareGlobalFunction("IsStdTableau");
DeclareGlobalFunction("CoxeterLength");
DeclareGlobalFunction("ContentVector");

DeclareGlobalFunction("SetDebugOutput");

# Some functionality for nice output:
DeclareGlobalFunction("CreateHTMLSummary");
DeclareGlobalFunction("PrintDecompositionMatrixAsLatex");
DeclareGlobalFunction("PrintBasisOfOrderAsMarkdown");
# And JupyterKernel integration, if available
if IsBound(JupyterRender) then
  DeclareGlobalFunction("JupyterDisplayDecompositionMatrix");
  DeclareGlobalFunction("JupyterDisplayBasisOfOrder");
  DeclareGlobalFunction("JupyterDisplayMultiMatrix");
fi;


# Internal use only
DeclareGlobalFunction("_QuickIso@");
DeclareGlobalFunction("_LatticeInvariants@");
DeclareGlobalFunction("_InverseRatMat@");
DeclareGlobalFunction("_Debug@");
DeclareGlobalFunction("_DebugWriteByte@");
DeclareGlobalVariable("_DEBUG_OUTPUT@");
