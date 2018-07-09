#############################################################################``
##  orders package - "orders.gi"
##  Copyright 2018 Florian Eisele
#############################################################################

# Types etc.

InstallGlobalFunction(_InverseRatMat@, function(m) return m^(-1); end);

InstallGlobalFunction(ZpOrderByMatrices, [IsPrime, IsList],
	function(p, gens)
 		local order;
		order := rec( p := p, gens := gens, n := Size(gens[1]) );
		return Objectify(NewType(ZpOrderFamily, IsZpOrder and IsZpOrderMatrixRep) , order);
  end
);

InstallGlobalFunction(ZpOrderByMultiMatrices, [IsPrime, IsList],
	function(p, gens)
		local order;
		order := rec( p := p, gens := gens, nvec := List(gens[1], Size) );
		return Objectify(NewType(ZpOrderFamily, IsZpOrder and IsZpOrderMultiMatrixRep) , order);
	end
);

InstallMethod(String, [IsZpOrder and IsZpOrderMatrixRep],
	function(order)
		return Concatenation("<order over Z", String(order!.p), " on ", String(Size(order!.gens)),
		      " generators in Z", String(order!.p), "^", String(order!.n), "x", String(order!.n),
		      ">");
	end
);

InstallMethod(String, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(order)
		return Concatenation("<order over Z", String(order!.p), " on ", String(Size(order!.gens)),
		      " generators in ", String(Size(order!.nvec)), " Wedderburn components>");
	end
);

InstallMethod(ViewString, [IsZpOrder and IsZpOrderMultiMatrixRep], String);

if IsBound(JupyterRenderable) then
	InstallMethod(JupyterRender, [IsZpOrder and IsZpOrderMultiMatrixRep],
		function(order)
			local str;
			str := Concatenation("&lt;order over $\\mathbb Z_{", String(order!.p), "}$ on ", String(Size(order!.gens)),
			      " generators in ", String(Size(order!.nvec)), " Wedderburn components&gt;");
			return JupyterRenderable(rec(("text/markdown") := str), rec());
		end
	);
fi;

InstallGlobalFunction(RModuleOverZpOrder, [IsZpOrder, IsList],
	function(order, rep)
		local modulus, n, M;

		# Check arguments:
		if Size(rep) <> Size(order!.gens) then
			Error("Wrong number of representation matrices");
		fi;
		if Size(Filtered(rep, m -> not IsMatrix(m))) <> 0 then
			Error("Second argument is supposed to be a list of matrices");
		fi;
		if rep[1][1][1] in Integers then
			modulus := 0;
			if Size(Filtered(rep, m -> not m in MatrixAlgebra(Integers, Size(rep[1])))) <> 0 then
				Error("Not all representation matrices are square matrices over the same ring");
			fi;
		elif rep[1][1][1] in GF(order!.p) then
			modulus := order!.p;
			if Size(Filtered(rep, m -> not m in MatrixAlgebra(GF(modulus), Size(rep[1])))) <> 0 then
				Error("Not all representation matrices are square matrices over the same ring");
			fi;
		else
			Error("Unsupported ground ring");
		fi;

		# Create the modules:
		n := Size(rep[1]);
		if modulus = 0 then
			M := Objectify(NewType(RModuleOverZpOrderFamily,
			                 IsRLatticeOverZpOrder and IsRModuleByRepresentationRep),
			                 rec( p := order!.p, n := n, gens := rep, order := order ) );
		else # modulus = order!.p
			M := Objectify(NewType(RModuleOverZpOrderFamily,
			                 IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep),
			                 rec( p := order!.p, n := n, gens := rep, order := order ) );
		fi;
		Setter(IsZero)(M, false);
		return M;
	end
);

InstallGlobalFunction(ZeroRModule, [IsZpOrder],
	function(order)
		local M;
		M := Objectify(NewType(RModuleOverZpOrderFamily,
	                IsRLatticeOverZpOrder and IsRModuleOverZpOrderModp and IsZeroRModuleRep),
	                rec( p := order!.p, order := order ) );
		Setter(IsZero)(M, true);
		return M;
	end
);

InstallMethod(ReduceModP, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(m)
		return RModuleOverZpOrder(m!.order, List(m!.gens, v -> One(GF(m!.order!.p))*v));
	end
);

InstallMethod(ReduceModP, [IsRModuleOverZpOrderModp],
	function(m)
		return m;
	end
);

InstallMethod(Dimension, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(m)
		return m!.n;
	end
);

InstallMethod(Dimension, [IsRModuleOverZpOrder and IsZeroRModuleRep],
	function(m)
		return 0;
	end
);

InstallMethod(Generators, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		return A!.gens;
	end
);

InstallMethod(Rep, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M)
		return M!.gens;
	end
);

InstallMethod(SimpleNames, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		if IsBound(A!.simple_names) then
			return A!.simple_names;
		else
			Error("The simple modules of this order were not named!");
		fi;
	end
);

InstallMethod(ComponentNames, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		if IsBound(A!.component_names) then
			return A!.component_names;
		else
			Error("The Wedderburn-components of this order were not named!");
		fi;
	end
);

InstallMethod(NameSimpleModules, [IsZpOrder and IsZpOrderMultiMatrixRep, IsList],
	function(A, v)
		if Size(v) <> Size(SimpleModules(A)) then
			Error("The second argument does not have the right  size!");
		elif not ForAll(v, IsString) then
			Error("The second argument is supposed to be a list of strings!");
		fi;
		A!.simple_names := StructuralCopy(v);
	end
);

InstallMethod(NameWedderburnComponents, [IsZpOrder and IsZpOrderMultiMatrixRep, IsList],
	function(A, v)
		if Size(v) <> Size(A!.nvec) then
			Error("The second argument does not have the right  size!");
		elif not ForAll(v, IsString) then
			Error("The second argument is supposed to be a list of strings!");
		fi;
		A!.component_names := StructuralCopy(v);
	end
);

InstallMethod(InstallTrivialEndomorphismRings, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		A!.Qp_end_bases := List(A!.nvec, n -> [IdentityMat(n, Integers)]);
		A!.Qp_end_bases_smallmats := List(A!.nvec, v -> [IdentityMat(1, Integers)]);
	end
);

InstallMethod(InstallEndomorphismRingsNC, [IsZpOrder and IsZpOrderMultiMatrixRep, IsList],
	function(A, E)
		local i, j;
		A!.Qp_end_bases := E;
		A!.Qp_end_bases_smallmats := List(E, x -> [ ]);
		# Calculate regular representations:
		# (TODO: Test this in a more non-trivial setting!)
		for i in [1..Size(E)] do
			for j in [1..Size(E[i])] do
				A!.Qp_end_bases_smallmats[i][j] := List(E[i], e -> Flat(e*E[i][j]))*
					TransposedMat(List(IdentityMat(Size(E[i])), v -> SolutionIntMat(TransposedMat(List(E[i], Flat)),
					v)));
			od;
		od;
	end
);

MakeReadWriteGlobal("_MAGMA_EXECUTABLE@");
_MAGMA_EXECUTABLE@ := "magma";

InstallGlobalFunction(SetMAGMAExecutable, [IsString],
	function(s)
		_MAGMA_EXECUTABLE@ := s;
	end
);

InstallMethod(CalculateEndomorphismRingsByReynoldsNC, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local endo;
		endo := List(IrreducibleLattices(A), L -> HomByReynoldsNC(L, L));
		InstallEndomorphismRingsNC(A, endo);
	end
);

InstallMethod(CalculateEndomorphismRingsWithMAGMA, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local fdin, fdout, out;
		fdin := TmpName();
		fdout := TmpName();
		out := OutputTextFile(fdin, false);
		SetPrintFormattingStatus(out, false);
		PrintTo(out, "m := ", A!.gens, ";\nret := [ ];\n");
		PrintTo(out, "PrintFile(\"", fdout, "\", \"local A;\\n A := \");\n");
		PrintTo(out, "for x in [1..", Size(A!.nvec),
			"] do\n\tV := RModule([ChangeRing(Matrix(v[x]), Integers()) : v in m]);\n",
			"\tAppend(~ret, AHom(V, V));\nend for;\nPrintFile(\"", fdout,
			"\", [ [RowSequence(qq) : qq in Basis(q)] : q in ret ]);\n");
		PrintTo(out, "PrintFile(\"", fdout, "\", \";\\nreturn A;\\n\");\n");
		PrintTo(out, "quit;\n");
		CloseStream(out);
		Exec(_MAGMA_EXECUTABLE@, "-b", fdin);
		Exec("rm", fdin);
		InstallEndomorphismRingsNC(A, ReadAsFunction(fdout)());
		Exec("rm", fdout);
	end
);

InstallMethod(String, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(m)
		return Concatenation("<lattice over an order given by representation on Z", String(m!.p),
		                   "^1x", String(m!.n),  ">");
	end
);

InstallMethod(ViewString, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep], String);

InstallMethod(String, [IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep],
	function(m)
		return Concatenation("<torsion module over an order given by representation on GF(",
		                    String(m!.p), ")^1x", String(m!.n),  ">");
	end
);

InstallMethod(ViewString, [IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep], String);

InstallMethod(String, [IsRModuleOverZpOrder and IsZeroRModuleRep],
	function(m)
		return Concatenation("<zero module over an order over Z", String(m!.p), ">");
	end
);

InstallMethod(ViewString, [IsRModuleOverZpOrder and IsZeroRModuleRep], String);


#InstallMethod(ViewObj, [IsList], 10,
#	function(l)
#		if Size(l) <= 10 then
#			TryNextMethod();
#		else
#			Print("<A list of size ", Size(l), ">");
#		fi;
#	end
#);

#InstallMethod(ViewObj, [IsMatrix],
#	function(l)
#		if Size(l) <= 10 and Size(l[1]) <= 10 then
#			TryNextMethod();
#		elif l[1][1] in Integers then
#			Print("<A ", Size(l), "x", Size(l[1]), " matrix over Z>");
#		else
#			Print("<A ", Size(l), "x", Size(l[1]), " matrix over ", Ring(l[1][1]), ">");
#		fi;
#	end
#);

InstallMethod(PrintAsFunction, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local comp, name, i;
		comp := Filtered(NamesOfComponents(A), v -> v <> "simple");
		Print("local A, B, S, i;\n");
		Print("A := rec( );\n");

		for name in comp do
			Print("A.", name, " := ", A!.(name), ";\n");
		od;
		Print("B := Objectify(NewType(ZpOrderFamily, IsZpOrder and IsZpOrderMultiMatrixRep), A);\n");

		if IsBound(A!.simple) then
			Print("S := [ ];");
			for i in [1..Size(A!.simple)] do
				Print("S[", i, "] := function()\n");
				PrintAsFunction(A!.simple[i]);
				Print("end;\n");
			od;
			Print("B!.simple := List(S, s -> s());\n");
			Print("for i in [1..Size(S)] do\n   B!.simple[i]!.order := B;\nod;\n");
		fi;

		Print("return B;\n");
	end
);

InstallMethod(PrintAsFunction, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M)
		local comp, name;
		comp := Filtered(NamesOfComponents(M), v -> v <> "order");
		Print("local M;\n");
		Print("M := rec( );\n");
		Print("M.order := \"NEEDS TO BE SET TO A VALUE!\";\n");
		for name in comp do
			Print("M.", name, " := ", M!.(name), ";\n");
		od;
		if IsRLatticeOverZpOrder(M) then
			Print("return Objectify(NewType(RModuleOverZpOrderFamily, IsRLatticeOverZpOrder and IsRModuleByRepresentationRep), M);\n");
		elif IsRModuleOverZpOrderModp(M) then
			Print("return Objectify(NewType(RModuleOverZpOrderFamily, IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep), M);\n");
		else
			Error("Unsupported module type! This should not happen!");
		fi;
	end
);

InstallMethod(SaveAsFunction, [IsString, IsObject],
	function(filename, X)
		PrintTo1(filename, function() PrintAsFunction(X); end);
	end
);

InstallMethod(SaveAsRecord, [IsString, IsZpOrder and IsZpOrderMultiMatrixRep],
	function(filename, A)
		local out;
		out := OutputTextFile(filename, false);
		SetPrintFormattingStatus(out, false);
		PrintTo(out, "local A;\nA := rec( );\n");
		PrintTo(out, "A.p := ", A!.p, ";\n");
		PrintTo(out, "A.nvec := ", A!.nvec, ";\n");
		PrintTo(out, "A.gens := ", A!.gens, ";\n");
		if IsBound(A!.component_names) then
			PrintTo(out, "A.component_names := ", A!.component_names, ";\n");
		fi;
		if IsBound(A!.simple) then
			PrintTo(out, "A.simple := ", List(A!.simple, v -> v!.gens), ";\n");
			PrintTo(out, "A.simple_multiplicities :=", A!.simple_multiplicities, ";\n");
		fi;
		if IsBound(A!.simple_names) then
			PrintTo(out, "A.simple_names := ", A!.component_names, ";\n");
		fi;
		if IsBound(A!.Qp_end_bases) then
			PrintTo(out, "A.Qp_end_bases := ", A!.Qp_end_bases, ";\n");
		fi;
		if IsBound(A!.idempot) then
			PrintTo(out, "A.idempot := ", A!.idempot, ";\n");
		fi;
		if IsBound(A!.blockmap) then
			PrintTo(out, "A.blockmap := ", A!.blockmap, ";\n");
		fi;
		if IsBound(A!.blockidx) then
			PrintTo(out, "A.blockidx := ", A!.blockidx, ";\n");
		fi;
		if IsBound(A!.gens_are_basis) then
			PrintTo(out, "A.gens_are_basis := ", A!.gens_are_basis, ";\n");
		fi;
		PrintTo(out, "return A;\n");
		CloseStream(out);
	end
);

InstallMethod(ReadOrder, [IsString],
	function(filename)
		return ReadAsFunction(filename)();
	end
);

InstallMethod(ReadModule, [IsString, IsZpOrder],
	function(filename, A)
		local M;
		M := ReadAsFunction(filename)();
		M!.order := A;
		return M;
	end
);

InstallGlobalFunction(UnFlattenMultiMatrixNC, [IsList, IsList],
	function(v, dimvec)
		local pos, v0, k;
		pos := List([1..Size(dimvec)], i -> Sum(List([1..i-1], j -> dimvec[j]^2)) + 1);
		v0 := List([1..Size(dimvec)], i -> v{[pos[i]..pos[i] + dimvec[i]^2 - 1]});
		for k in [1..Size(dimvec)] do
			v0[k] := List([1..dimvec[k]], i -> v0[k]{[(i-1)*dimvec[k]+1..i*dimvec[k]]});
		od;
		return v0;
	end
);

InstallMethod(DirectSumOfModules, [IsList],
	function(l)
		local m;
		if Size(l) = 0 then
			Error("Expecting a list of at least one module");
		fi;
		if ForAll(l, IsRLatticeOverZpOrder) or ForAll(l, IsRModuleOverZpOrderModp) then
			if ForAny(l, m -> m!.order <> l[1]!.order) then
				Error("Modules must all be defined over the same order!");
			fi;
			m := RModuleOverZpOrder(l[1]!.order, List([1..Size(l[1]!.gens)], n -> DiagonalJoin(List(l, l -> l!.gens[n]))));
			if ForAll(l, x -> IsBound(x!.embedding_into_irr_lat)) then
				m!.embedding_into_irr_lat := [DiagonalJoin(List(l, x -> x!.embedding_into_irr_lat[1])),
				                              Concatenation(List(l, x -> x!.embedding_into_irr_lat[2]))];
			fi;
			return m;
		else
			TryNextMethod();
		fi;
	end
);

InstallMethod(DirectSumOfModules, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep,
                                   IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M, N)
		return DirectSumOfModules([M, N]);
	end
);

InstallMethod(IrreducibleLattices, [IsZpOrder and IsZpOrderMatrixRep],
	function(A)
		local M;
		M := RModuleOverZpOrder(A, A!.gens);
		M!.embedding_into_irr_lat := [IdentityMat(M!.n, Integers), [1]];
		if IsBound(A!.Qp_end_basis) then
			M!.Qp_end_basis := A!.Qp_end_basis;
		fi;
		if IsBound(A!.Qp_end_basis_smallmats) then
			M!.Qp_end_basis_smallmats := A!.Qp_end_basis_smallmats;
		fi;
		return [M];
	end
);

InstallMethod(IrreducibleLattices, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local M, i;
		M := List([1..Size(A!.nvec)], i -> RModuleOverZpOrder(A, List(A!.gens, g -> g[i])));
		for i in [1..Size(M)] do
			M[i]!.embedding_into_irr_lat := [IdentityMat(M[i]!.n, Integers), [i]];
			if IsBound(A!.Qp_end_bases) then
				if IsBound(A!.Qp_end_bases[i]) then
					M[i]!.Qp_end_basis := A!.Qp_end_bases[i];
				fi;
			fi;
			if IsBound(A!.Qp_end_bases_smallmats) then
				if IsBound(A!.Qp_end_bases_smallmats[i]) then
					M[i]!.Qp_end_basis_smallmats := A!.Qp_end_bases_smallmats[i];
				fi;
			fi;
		od;
		return M;
	end
);

InstallMethod(IrreducibleLattices, [IsZpOrder and IsZpOrderMatrixRep],
	function(A)
		local M;
		M := RModuleOverZpOrder(A, A!.gens);
		M!.embedding_into_irr_lat := [IdentityMat(M!.n, Integers), [1]];
		if IsBound(A!.Qp_end_basis) then
			M!.Qp_end_basis := A!.Qp_end_basis;
		fi;
		if IsBound(A!.Qp_end_basis_smallmats) then
			M!.Qp_end_basis_smallmats := A!.Qp_end_basis_smallmats;
		fi;
		return [M];
	end
);

InstallMethod(BlocksOfZpOrder, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local M, blocks, curr, last, newlast, ahead;
		M := DecompositionMatrix(A);
		ahead := [1..Size(M)];
		last := [ ];
		blocks := [ ];
		while not IsZero(Size(ahead)) do
			if IsZero(Size(last)) then
				last := [ Remove(ahead, 1) ];
				curr := [ ];
			else
				last :=  Filtered(ahead, i -> not IsZero(List(last, j -> M[j])*TransposedMat([M[i]])));
				ahead := Filtered(ahead, i -> not i in last);
			fi;
			Append(curr, last);
			if IsZero(Size(last)) or IsZero(Size(ahead)) then
				Sort(curr);
				Add(blocks,curr);
			fi;
		od;
		return blocks;
	end
);

InstallMethod(ExtractWedderburnComponents, [IsZpOrder and IsZpOrderMultiMatrixRep, IsList],
	function(A, l)
		local B, s, i;
		if not IsDuplicateFree(l) then
			Error("Second argument is supposed to be a duplicate-free list!");
		elif not IsSubset(AsSet([1..Size(A!.nvec)]), AsSet(l)) then
			Error("Entries of the second argument are supposed to be indices of Wedderburn-components!");
		fi;

		B := ZpOrderByMultiMatrices(A!.p, List(A!.gens, g -> g{l}));
		if IsBound(A!.component_names) then
			B!.component_names := A!.component_names{l};
		fi;
		if IsBound(A!.Qp_end_bases) then
			B!.Qp_end_bases := A!.Qp_end_bases{l};
		fi;
		if IsBound(A!.Qp_end_bases_smallmats) then
			B!.Qp_end_bases_smallmats := A!.Qp_end_bases_smallmats{l};
		fi;
		if IsBound(A!.simple) and IsBound(A!.simple_multiplicities) then
			s := Filtered([1..Size(A!.simple)], i -> not IsZero(A!.simple_multiplicities{l}{[i]}));
			B!.simple := StructuralCopy(A!.simple{s});
			for i in [1..Size(B!.simple)] do
				B!.simple[i]!.order := B;
			od;
			B!.simple_multiplicities := A!.simple_multiplicities{l}{s};
			if IsBound(A!.simple_names) then
				B!.simple_names := A!.simple_names{s};
			fi;
			if IsBound(A!.idempot) then
				B!.idempot := A!.idempot{s};
			fi;
			if IsBound(A!.blockmap) then;
				B!.blockmap := Filtered(A!.blockmap, x -> x[1][1] in s and x[1][2] in s);
			fi;
			if IsBound(A!.blockidx) then
				B!.blockidx := A!.blockidx;
			fi;
			B!.gens_are_basis := false;
		fi;
		if IsBound(A!.eval_map) then
			B!.eval_map := A!.eval_map{l};
		fi;
		return B;
	end
);

InstallMethod(DirectSumOfOrders, [IsList],
	function(l)
		local B, Bgens, g, zero, i, j;
		if IsZero(Size(l)) then
			Error("Argument list empty!");
		elif not ForAll(l, IsZpOrderMultiMatrixRep) then
			Error("Expects a list of orders!");
		elif ForAny(l, z -> z!.p <> l[1]!.p) then
			Error("Not all orders defined in same characteristic!");
		fi;

		zero := List(l, z -> List(z!.nvec, n -> NullMat(n, n)));
		Bgens := [ ];
		for i in [1..Size(l)] do
			g := zero;
			for j in [1..Size(l[i]!.gens)] do
				g[i] := l[i]!.gens[j];
				Add(Bgens, Concatenation(g));
			od;
		od;
		B := ZpOrderByMultiMatrices(l[1]!.p, Bgens);

		if ForAll(l, z -> IsBound(z!.gens_are_basis)) then
			B!.gens_are_basis := ForAll(l, z -> z!.gens_are_basis);
		fi;
		if ForAll(l, z -> IsBound(z!.Qp_end_bases)) then
			B!.Qp_end_bases := Concatenation(List(l, z -> z!.Qp_end_bases));
		fi;
		if ForAll(l, z -> IsBound(z!.Qp_end_bases_smallmats)) then
			B!.Qp_end_bases_smallmats := Concatenation(List(l, z -> z!.Qp_end_bases_smallmats));
		fi;
		if ForAll(l, z -> IsBound(z!.component_names)) then
			B!.component_names := Concatenation(List(l, z -> z!.component_names));
		fi;
		# [...TODO... (maybe) transfer simple modules]

		return B;
	end
);

InstallMethod(DirectSumOfOrders, [IsZpOrder, IsZpOrder],
	function(A, B)
		return DirectSumOfOrders([A, B]);
	end
);

InstallMethod(NameSimpleModulesByDims, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local dims, names;
		dims := List(SimpleModules(A), s -> s!.n);
		names := List([1..Size(dims)], i -> Concatenation(String(dims[i]), Concatenation(List(Filtered(dims{[1..i-1]},
			j -> j = dims[i]), z -> "'"))));
		A!.simple_names := names;
	end
);

InstallMethod(NameWedderburnComponentsByDims, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local dims, names;
		dims := A!.nvec;
		names := List([1..Size(dims)], i -> Concatenation(String(dims[i]), Concatenation(List(Filtered(dims{[1..i-1]},
			j -> j = dims[i]), z -> "'"))));
		A!.component_names := names;
	end
);

## Algorithms

InstallMethod(SimpleModules, [IsZpOrder and IsZpOrderMatrixRep],
	function(order)
		local F, M, simple, col;

		if IsBound(order!.simple) then
			return order!.simple;
		fi;

		F := GF(order!.p);
		M := GModuleByMats(List(order!.gens, m -> One(F)*m), F);
		col := MTX.CollectedFactors(M);
		simple := List(col, m -> RModuleOverZpOrder(order,MTX.Generators(m[1])));
		order!.simple_multiplicities := List(col, m -> m[2]);

		order!.simple := simple;
		return simple;
	end
);

InstallMethod(SimpleModules, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(order)
		local F, M, simple, mult, col, i, j, k, isnew;

		if IsBound(order!.simple) then
			return order!.simple;
		fi;

		F := GF(order!.p);
		M := List([1..Size(order!.nvec)], z -> GModuleByMats(List(order!.gens, m -> One(F)*m[z]), F));
		col := List(M, MTX.CollectedFactors);

		simple := List(col[1], z -> z[1]);
		mult := [List(col[1], z -> z[2])];
		for i in [2..Size(col)] do
			Add(mult, [ ]);
			for j in [1..Size(col[i])] do
				isnew := true;
				for k in [1..Size(simple)] do
					if Size(MTX.Homomorphisms(simple[k], col[i][j][1])) <> 0 then
						isnew := false;
						mult[i][k] := col[i][j][2];
						break;
					fi;
				od;
				if isnew then
					Add(simple, col[i][j][1]);
					mult[i][Size(simple)] := col[i][j][2];
				fi;
			od;
		od;

		for i in [1..Size(mult)] do
			for j in [1..Maximum(List(mult, Size))] do
				if not IsBound(mult[i][j]) then
					mult[i][j] := 0;
				fi;
			od;
		od;

		simple := List(simple, m -> RModuleOverZpOrder(order, MTX.Generators(m)));
		order!.simple := simple;
		order!.simple_multiplicities := mult;

		return simple;
	end
);

InstallMethod(DecompositionMatrix, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		if not IsBound(A!.simple_multiplicities) then
			SimpleModules(A);
		fi;
		return A!.simple_multiplicities;
	end
);

InstallMethod(DecompositionMatrix, [IsZpOrder and IsZpOrderMatrixRep],
	function(A)
		if not IsBound(A!.simple_multiplicities) then
			SimpleModules(A);
		fi;
		return [ A!.simple_multiplicities ];
	end
);

InstallMethod(MaximalSubmoduleBases, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep,
                                      IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep],
	function(M, S) # Determines maximal submodules as kernels of homomorphisms onto S
		local hom0, hom1, hom2, im, k, V, im2, max, m, E;

		hom0 := HomToSimpleNC(M, S);
		if Size(hom0) = 0 then
			return [ ];
		fi;
		im := VectorSpace(GF(M!.p), TransposedMat(hom0[1]));
		hom1  := [hom0[1]];
		for k in [2..Size(hom0)] do
			V := VectorSpace(GF(M!.p), TransposedMat(hom0[k]));
			im2 := im + V;
			if Dimension(im2) > Dimension(im) then
				Add(hom1, hom0[k]);
				im := im2;
			fi;
		od;
		Unbind(hom0);

		# hom1 now is an End(S)-vectorspace-basis of Hom(M, S). Now we calculate
		# representatives for all End(S)-subspaces of Hom(M, S)
		E := VectorSpace(GF(M!.p), HomToSimpleNC(S, S), "basis");
		hom2 := [ hom1[1] ];
		for k in [2..Size(hom1)] do
			hom2 := Concatenation(Concatenation(List(E, e -> List(hom2, h -> h + hom1[k]*e))), [hom1[k]]);
		od;
		Unbind(hom1);
		max := List(hom2, h -> NullspaceMat(h));

		# max now contains bases of the kernels over GF(p). Now we have to
		# construct the corresponding modules:
		if IsRLatticeOverZpOrder(M) then
			for m in [1..Size(max)] do
				max[m] := List(max[m], x -> IntVecFFE(x));
				Append(max[m], M!.p*IdentityMat(M!.n, Integers));
				max[m] := BaseIntMat(max[m]);
			od;
		elif IsRModuleOverZpOrderModp(M) then
			# We have to filter out the zero-modules:
			# max := Filtered(max, m -> Size(m) <> 0);
		else
			Error("Unexpected (and unsupported) module type. This shouldn't happen.");
		fi;

		return max;
	end
);

InstallMethod(MaximalSubmoduleBases, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M) # Calculate all maximal submodules as kernels of homomorphisms onto simple modules
               # and return tuples [N, S], with M < N maximal and M/N=S
		local simple, S, max;

		simple := SimpleModules(M!.order);
		max := [ ];
		for S in simple do
			Append(max, List(MaximalSubmoduleBases(M, S), m -> [m, S]));
		od;
		return max;
	end
);

InstallMethod(SubmoduleByBasisNC, [IsRModuleOverZpOrder, IsList],
	function(M, b)
		if IsZero(Size(b)) then
			return ZeroRModule(M!.order);
		else
			TryNextMethod();
		fi;
	end
);

InstallMethod(SubmoduleByBasisNC, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep, IsMatrix],
	function(M, b)
			local m, binv, emb2, V, W, j, pos, d, idx, idx2;

			if IsRLatticeOverZpOrder(M) and (Size(b) = M!.n) then
				# In the full sublattice case we just conjugate the representation
				binv := _InverseRatMat@(b);
				m := RModuleOverZpOrder(M!.order, List(M!.gens, x -> b*x*binv));
				if IsBound(M!.embedding_into_irr_lat) then
					m!.embedding_into_irr_lat := [b*M!.embedding_into_irr_lat[1], M!.embedding_into_irr_lat[2]];
				fi;
			elif IsRLatticeOverZpOrder(M) then
				binv := RightInverse(b);
				m := RModuleOverZpOrder(M!.order, List(M!.gens, x -> b*x*binv));
				if IsBound(M!.embedding_into_irr_lat) then # This part is essentially UNTESTED!
					emb2 := b*M!.embedding_into_irr_lat[1]; # This needs to be composed with an adequate projection
					V := Subspace(FullRowSpace(Rationals, Size(b)), [ ]);
					pos := 1;
					idx := [ ];
					idx2 := [ ];
					for j in M!.embedding_into_irr_lat[2] do
						d := M!.order!.nvec[j];
						W := Subspace(FullRowSpace(Rationals, Size(b)), TransposedMat(emb2){[pos..pos+d-1]});
						if not IsSubspace(V, W) then
							Append(idx, [pos..pos+d-1]);
							Add(idx2, j);
							V := V + W;
						fi;
						pos := pos + d;
					od;
					m!.embedding_into_irr_lat := [ List(emb2, v -> v{idx}), idx2 ];
				fi;
			elif IsRModuleOverZpOrderModp(M) then
				# In the torsion module case we need a right inverse of the basis matrix,
				# since it doesn't have to be square
				binv := RightInverse(b);
				m := RModuleOverZpOrder(M!.order, List(M!.gens, x -> b*x*binv));
			else
				Error("Unexpected (and unsupported) module type. This shouldn't happen.");
			fi;
			return m;
	end
);

InstallMethod(MaximalSubmodules, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep,
                                  IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep],
	function(M, S) # Returns a list of maximal submodules of M with quotient S
		local bases, b, max;

		bases := MaximalSubmoduleBases(M, S);
		max := [  ];
		for b in bases do
			Add(max, SubmoduleByBasisNC(M, b));
		od;

		return max;
	end
);

InstallMethod(MaximalSubmodules, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M) # Returns pairs [N, S], where N is a maximal submodule of M and S = M/N
		local simple, S, max;

		simple := SimpleModules(M!.order);
		max := [ ];
		for S in simple do
			Append(max, List(MaximalSubmodules(M,S), m -> [m, S]));
		od;

		return max;
	end
);

InstallMethod(MaximalSubmoduleBasesMTX, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M) # Returns a list [B_1, ... B_k] of maximal submodules given by basis-matrices
		local MM, max, m;

		MM := GModuleByMats(List(M!.gens, g -> One(GF(M!.p))*g), GF(M!.p));
		max := MTX.BasesMaximalSubmodules(MM);

		if IsRModuleOverZpOrderModp(M) then
			return max;
		elif IsRLatticeOverZpOrder(M) then
			for m in [1..Size(max)] do
				max[m] := List(max[m], x -> IntVecFFE(x));
				Append(max[m], M!.p*IdentityMat(M!.n, Integers));
				max[m] := BaseIntMat(max[m]);
			od;
			return max;
		else
			Error("Unexpected (and unsupported) module type. This shouldn't happen.");
		fi;
	end
);

InstallGlobalFunction(LatticeAlgorithm,
	[IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(M, b) # b = false => Find representatives for all iso.classes of lattices in M
                  # b = true => Just find the ones with simple radical quotient
		local simple, lat, lastlayer, thislayer, m, n, x, max, bas;

		simple := SimpleModules(M!.order);
		lat := [ ];

		lastlayer := [ IdentityMat(M!.n, Integers) ];

		while Size(lastlayer) > 0 do
			_Debug@("New layer, found so far: ", Size(lat), "; looking at ", Size(lastlayer), " lattices", "\n");
			thislayer := [ ];
			for m in lastlayer do
				max := MaximalSubmoduleBases(SubmoduleByBasisNC(M, m));

			   # Add m to the list that will be returned (if the conditions specified by b are met)
				if ((Size(max) = 1) and b) then
					Append(lat, [[m, max[1][2]]]);
				elif not b then
					Append(lat, [m]);
				fi;

				# Test which lattices are new and not contained in p*M
				for n in max do
					bas := n[1]*m;
					if  IsZero(bas mod M!.p) then # Lattice is contained in p*M
						continue;
					fi;
					NormalFormIntMat(bas, 0+2+16);
					if ForAny(thislayer, x -> x = bas) then
						continue;
					fi;
					Add(thislayer, bas);
				od;
			od;
			lastlayer := thislayer;
		od;

		return lat;
	end
);

InstallGlobalFunction(LatticesWithSimpleRadQuo, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(M)
		local lat;

		lat := LatticeAlgorithm(M, true);
		return List(lat, l -> [SubmoduleByBasisNC(M, l[1]), l[2]]);
	end
);

InstallGlobalFunction(GlueUpNC, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep,
                                 IsRLatticeOverZpOrder and IsRModuleByRepresentationRep,
											IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep],
	function(P, Q, S)
		local alpha, beta, betal, H, JQ, T, X, cont, L, max, N, i, j;

		alpha := HomToSimpleNC(P, S)[1];
		beta  := HomToSimpleNC(Q, S)[1];
		betal := List(IdentityMat(S!.n, GF(S!.p)), e -> SolutionMat(beta, e));
		H := List(alpha*betal, z -> List(z, zz -> Int(zz)));
		JQ := MaximalSubmoduleBases(Q, S)[1];
		T := MutableNullMat(P!.n + Q!.n, P!.n + Q!.n, Integers);
		for i in [1..P!.n] do
			T[i][i] := 1;
			for j in [1..Q!.n] do
				T[i][P!.n + j] := H[i][j];
			od;
		od;
		for i in [1..Q!.n] do
			for j in [1..Q!.n] do
				T[P!.n + i][P!.n + j] := JQ[i][j];
			od;
		od;

		X := RModuleOverZpOrder(P!.order, List([1..Size(P!.gens)], n -> DiagonalJoin([P!.gens[n], Q!.gens[n]])));
		X := SubmoduleByBasisNC(X, T);

		repeat
			cont := false;
			for L in SimpleModules(S!.order) do
				max := MaximalSubmoduleBases(X, L);
				for N in max do
					if Rank(One(GF(S!.p)) * N{[1..Size(N)]}{[1..P!.n]}) = P!.n then
						cont := true;
						T := N*T;
						X := SubmoduleByBasisNC(X, N);
						break;
					fi;
				od;
				if cont then
					break;
				fi;
			od;
			_DebugWriteByte@(STDOut, 43);
		until not cont;
		_Debug@("\n");

		if IsBound(P!.embedding_into_irr_lat) and IsBound(Q!.embedding_into_irr_lat) then
			X!.embedding_into_irr_lat :=
			  [T*DiagonalJoin(P!.embedding_into_irr_lat[1], Q!.embedding_into_irr_lat[1]),
				Concatenation(P!.embedding_into_irr_lat[2], Q!.embedding_into_irr_lat[2])];
		fi;
		return [X, T];
	end
);

InstallMethod(ProjectiveIndecomposableLattices, [IsZpOrder and IsZpOrderMatrixRep],
	function(A)
		local M, S, P0, Pk, k, d, i, X, Xb, XT, ret;

		if not IsBound(A!.Qp_end_basis) then
			Error("Endomorphism rings of irreducible representation must be known!");
		fi;

		M := IrreducibleLattices(A)[1];
		S := SimpleModules(A);
		ret := [ ];
		P0 := LatticeAlgorithm(M, true);
		for k in [1..Size(S)] do
			d := A!.simple_multiplicities[k] * Size(HomToSimpleNC(S[k], S[k])) / Size(A!.Qp_end_basis);
                    # Multiplicity of Qp(x)M in the projective cover of S[k] (according to Brauer reciprocity)
			# if d > 2 or Size(A!.Qp_end_basis) > 1 then
			# 	Print("# Warning: Decomposition number > 2. This is not guaranteed to terminate!\n");
			# fi;
			Pk := List(Filtered(P0, m -> Size(HomToSimpleNC(m[2], S[k])) <> 0), l -> l[1]);
			X := SubmoduleByBasisNC(M, Pk[1]);
			for i in [2..d] do
				XT := GlueUpNC(SubmoduleByBasisNC(M, Pk[i]), X, S[k]);
				X := XT[1];
			od;
			Add(ret, [X, S[k]]);
			_Debug@("Done with ", S[k], "\n");
		od;

		return ret;
	end
);

InstallMethod(ProjectiveIndecomposableLattices, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local W, w, D, S, s, i, P0, P0s, P0w, P, Ps, m;
		if IsBound(A!.pim) then
			return List([1..Size(SimpleModules(A))], x -> [A!.pim[x], A!.simple[x]]);
		fi;
		S := SimpleModules(A);
		D := DecompositionMatrix(A);
		W := List([1..Size(A!.nvec)], n -> ZpOrderByMatrices(A!.p, List(A!.gens, x -> x[n])));
		for i in [1..Size(W)] do
			W[i]!.simple := A!.simple{Filtered([1..Size(A!.simple)], n -> D[i][n] <> 0)};
			W[i]!.simple_multiplicities := D[i]{Filtered([1..Size(A!.simple)], n -> D[i][n] <> 0)};
			W[i]!.Qp_end_basis := A!.Qp_end_bases[i];
		od;
		P0 := [ ];
		for i in [1..Size(W)] do
			P0w := ProjectiveIndecomposableLattices(W[i]);
			for m in P0w do
				m[1]!.order := A;
				m[1]!.embedding_into_irr_lat[2] := List(m[1]!.embedding_into_irr_lat[2], n -> i);
				m[2]!.order := A;
			od;
			Append(P0, P0w);
			_Debug@("***> Done with component\n");
		od;
		P := [ ];
		for s in S do
			P0s := List(Filtered(P0, x -> Size(HomToSimpleNC(s, x[2])) <> 0), y -> y[1]);
			Ps := P0s[1];
			for i in [2..Size(P0s)] do
				_Debug@("Glueing up: Step ", i-1, " in ", Size(P0s)-1, "\n");
				Ps := GlueUpNC(P0s[i], Ps, s)[1];
			od;
			Add(P, [Ps, s]);
			_Debug@("***> Done with ", s, "\n");
		od;

		A!.pim := List(P, x -> x[1]);
		return P;
	end
);

InstallMethod(BasicOrder, [IsList],
	function(PS) # Calculates (+)_i (+)_j Hom(PS[1][i],PS[1][j])^tr, where PS is a list with
	             # entries of the form [projective lattice, head of the aforementioned]
		local A, P, S, dimvec, homs, support, i, j, n, m, m0, posi, posj, bas, nf,
		      T, Tinv, cont, idempot, blockidx, blockmap, bas0, B, X, endo, endo_small, simple_, simple_mul_;
		if ForAny(PS, z -> Size(z) <> 2) or IsZero(Size(PS)) then
			Error("Expecting non-empty list of tuples [projective lattice, simple module]!");
		fi;
		P := List(PS, z -> z[1]);
		S := List(PS, z -> z[2]);
		if ForAny(P, p -> not IsRLatticeOverZpOrder(p)) then
			Error("First entries are supposed to be of lattices over an order!");
		fi;
		if ForAny(S, p -> not IsRModuleOverZpOrderModp(p)) then
			Error("Second entries are supposed to be torsion modules over an order!");
		fi;
		A := P[1]!.order;
		if ForAny(P, p -> p!.order <> A) or ForAny(S, p -> p!.order <> A) then
			Error("Not all modules are defined over the same order!");
		fi;

		homs := List(P, p -> List(P, q -> HomForLattices(p, q, [2])[1]));
		dimvec := List([1..Size(IrreducibleLattices(A))],
		               n -> List([1..Size(P)], function(i) if IsBound(homs[i][i][1][n]) then
		                                                   	return Size(homs[i][i][1][n]);
		                                                   else
		                                                      return 0; fi; end));
		support := Filtered([1..Size(dimvec)], z -> not IsZero(Sum(dimvec[z])));
		dimvec := dimvec{support};
		bas := [ ];
		idempot := [ ];
		blockidx := [ ];
		blockmap := [ ];
		for i in [1..Size(P)] do
			for j in [1..Size(P)] do
				if i = j then
					idempot[i] := List(dimvec, d -> DiagonalJoin(List(Filtered([1..Size(P)], zz->d[zz] <> 0),
							function(k) if k = i then
					          return IdentityMat(d[k], Integers);
					      else return NullMat(d[k], d[k], Integers); fi; end)));
				fi;
				if not IsZero(Size(homs[i][j])) then
					Add(blockmap, [[i, j], Size(blockidx) + 1]);
					Add(blockidx, [Size(bas) + 1, Size(homs[i][j])]);
				fi;
				for m in homs[i][j] do
					m0 := List([1..Size(dimvec)], k -> MutableNullMat(Sum(dimvec[k]), Sum(dimvec[k]), Integers));
					for n in [1..Size(dimvec)] do
						if IsBound(m[support[n]]) then
							posi := Sum(dimvec[n]{[1..i-1]});
							posj := Sum(dimvec[n]{[1..j-1]});
							m0[n]{[1..dimvec[n][i]] + posi}{[1..dimvec[n][j]] + posj} := m[support[n]];
						fi;
					od;
					Add(bas, List(m0, TransposedMat));
				od;
			od;
		od;

		endo_small := List(support, z -> A!.Qp_end_bases_smallmats[z]);
		endo := List([1..Size(support)], z -> List(endo_small[z], zz -> DiagonalJoin(List([1..Sum(dimvec[z])/Size(zz)], z_ -> zz))));

		# The matrices don't have to be integral at this point. Hence we apply an integral spinning algorithm:
		for i in [1..Size(dimvec)] do
			T := IdentityMat(Sum(dimvec[i]), Integers);
			repeat
				T := Concatenation(T, Concatenation(List(bas, b -> b[i])));
				T := BaseIntMat(Lcm(List(Flat(T), DenominatorRat)) * T);
				# Get the matrices in a (block-)diagonal form: (the ordering of the idempotents is important here)
				T := Concatenation(List(idempot, e -> BaseIntMat(T*e[i])));
				Tinv := T^-1;
				endo[i] := List(endo[i], z -> T*z*Tinv);
				cont := false;
				for j in [1..Size(bas)] do
					bas[j][i] := T*bas[j][i]*Tinv;
					cont := cont or not (bas[j][i] in Integers^[Sum(dimvec[i]), Sum(dimvec[i])]);
				od;
			until not cont;
		od;

		# Now we do some cosmetics (smaller entries):
		for i in blockidx do
			bas0 := List(bas{[i[1]..i[1]+i[2]-1]}, Flat);
			nf := NormalFormIntMat(bas0, 1+4+8);
			for j in [1..Size(nf.normal)] do
				if not IsZero(nf.normal[j][j]) then
					nf.normal[j][j] := A!.p^Valuation(nf.normal[j][j], A!.p);
				fi;
			od;
			bas0 := HermiteNormalFormIntegerMat(nf.rowtrans^-1 * nf.normal * nf.coltrans^-1);
			bas{[i[1]..i[1]+i[2]-1]} := List(bas0, b -> UnFlattenMultiMatrixNC(b, List(dimvec, Sum)));
		od;

		B := ZpOrderByMultiMatrices(A!.p, bas);
		B!.idempot := idempot;
		B!.blockmap := blockmap;
		B!.blockidx := blockidx;
		B!.Qp_end_bases := endo;
		B!.Qp_end_bases_smallmats := endo_small;
		if IsBound(A!.component_names) then
			B!.component_names := A!.component_names{support};
		fi;

		# Determine, order and name the simple B-modules;
		X := SimpleModules(B);
		simple_ := [ ];
		simple_mul_ := [ ];
		if IsBound(A!.simple_names) then
			B!.simple_names := [ ];
		fi;
		for j in [1..Size(S)] do
			for i in [1..Size(X)] do
				# Here we use that the first basis element in Hom(P[j],P[j])^tr in B!.gens
				# is not contained in the radical. That has been ensured by calculating Hermite-forms above.
				if not IsZero(X[i]!.gens[blockidx[Filtered(blockmap, z -> z[1] = [j, j])[1][2]][1]][1]) then
					simple_[j] := X[i];
					simple_mul_[j] := TransposedMat(B!.simple_multiplicities)[i];
					if IsBound(A!.simple_names) then
						B!.simple_names[j] := A!.simple_names[Filtered([1..Size(SimpleModules(A))], v ->
							Size(HomToSimpleNC(S[j], SimpleModules(A)[v])) <> 0)[1]];
					fi;
				fi;
			od;
		od;
		B!.simple := simple_;
		B!.simple_multiplicities := TransposedMat(simple_mul_);

		B!.gens_are_basis := true;

		return B;
	end
);

InstallMethod(BasicOrder, [IsZpOrder],
	function(A)
		return BasicOrder(ProjectiveIndecomposableLattices(A));
	end
);

InstallMethod(RadicalOfModule, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M)
		local S, hom2S, jS, max, edim, rad, J;
		S := SimpleModules(M!.order);
		hom2S := List(S, s -> HomToSimpleNC(M, s));
		max := List(Filtered(hom2S, v -> Size(v) <> 0), z -> NullspaceMat(TransposedMat(BaseMat(Concatenation(List(z, TransposedMat))))));
		edim := List(S, s -> Size(HomToSimpleNC(s,s)));
		rad := AsList(Basis(Intersection(List(max, v -> VectorSpace(GF(M!.p), v, Zero(FullRowSpace(GF(M!.p), M!.n)))))));
		if IsRLatticeOverZpOrder(M) then
			J := SubmoduleByBasisNC(M, BaseIntMat(Concatenation(List(rad, IntVecFFE), M!.p * IdentityMat(M!.n))));
		elif IsRModuleOverZpOrderModp(M) then
			J := SubmoduleByBasisNC(M, rad);
		else
			TryNextMethod();
		fi;
		return [J, List([1..Size(edim)], i -> Size(hom2S[i])/edim[i])];
	end
);

InstallMethod(RadicalOfModule, [IsRModuleOverZpOrder and IsZeroRModuleRep],
	function(M)
		return [M, List([1..Size(SimpleModules(M!.order))], v -> 0)];
	end
);

InstallMethod(RadicalSeries, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep, IsPosInt],
	function(M, n)
		local J, Xi, lst, i;
		J := M;
		lst := [ ];
		for i in [1..n] do
			Xi := RadicalOfModule(J);
			J := Xi[1];
			Add(lst, Xi[2]);
		od;
		return lst;
	end
);

InstallMethod(TopEpimorphism, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M) # Returns [M/Rad M,phi,v,hom], where phi is an epimorphism from M onto M/Rad M. v Is a vector of
			# multiplicities (i. e., v[i] is the multiplicity of the i-th simple module in M/Rad M).
			# hom is a list of End(S)-bases of Hom(M, S) (for all simple modules S)
		local i, k, V, im, hom0, hom1, R, phi;
		hom1 := [ ];
		for i in [1..Size(SimpleModules(M!.order))] do
			hom0 := HomToSimpleNC(M, SimpleModules(M!.order)[i]);
			if IsZero(Size(hom0)) then
				hom1[i] := [ ];
				continue;
			fi;
			im := VectorSpace(GF(M!.p), TransposedMat(hom0[1]));
			hom1[i]  := [hom0[1]];
			for k in [2..Size(hom0)] do
				V := VectorSpace(GF(M!.p), TransposedMat(hom0[k]));
				if not IsSubspace(V, im) then
					Add(hom1[i], hom0[k]);
					im := im + V;
				fi;
			od;
		od;
		# hom1[i] now is an End(S[i])-vectorspace-basis of Hom(M, S[i]]).

		R := DirectSumOfModules(Flat(List([1..Size(SimpleModules(M!.order))], l -> List(hom1[l],
			q -> SimpleModules(M!.order)[l]))));
		phi := List([1..M!.n], l -> Flat(List(hom1, q -> List(q, v -> v[l]))));
		return [R, phi, List(hom1, Size), hom1];
	end
);

InstallMethod(LiftHomomorphismNC, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep, IsMatrix,
	IsRModuleOverZpOrder and IsRModuleByRepresentationRep, IsMatrix,
	IsRModuleOverZpOrder],
	function(P, phi, M, psi, N) # This lifts (if possible) the homomorphism phi: P -> N to a homorphism
	                            # mu: P -> M such that mu.psi: P -> N is equal to phi.
		local hom0, hom, coeff, ret;
		if IsZero(N) then
			return MutableNullMat(P!.n, M!.n, DefaultFieldOfMatrix(M!.gens[1]));
		elif not IsRModuleByRepresentationRep(N) then
			TryNextMethod();
		fi;
		hom0 := Hom(P, M);
		hom := List(hom0, v -> Flat(v*psi));
		if DefaultFieldOfMatrix(hom) = Rationals then
			coeff := SolutionIntMat(hom, Flat(phi));
		else
			coeff := SolutionMat(hom, Flat(phi));
		fi;
		if coeff = fail then
			return fail;
		fi;
		if IsRModuleOverZpOrderModp(N) then
			coeff := IntVecFFE(coeff);
		fi;
		ret := Sum([1..Size(hom)], i -> coeff[i]*hom0[i]);
		return ret;
	end
);

InstallMethod(LiftHomomorphismToDirectSumNC, [IsList, IsList, IsMatrix,
	IsRModuleOverZpOrder and IsRModuleByRepresentationRep, IsMatrix,
	IsRModuleOverZpOrder],
	function(lst, svec, phi, M, psi, N) # This lifts (if possible) the homomorphism phi: P -> N to a homorphism
	                            # mu: P -> M such that mu.psi: P -> N is equal to phi.
	                            # (here P is defined as DirectSum(List(svec, lst[i])))
		local hom0, hom, coeff, ret, phom, i, h, hx, dimP, px, hx_;
		if not ForAll(lst, IsRModuleByRepresentationRep and IsRLatticeOverZpOrder) then
			TryNextMethod();
		fi;
		dimP := Sum(List(svec, q -> lst[q]!.n));
		if IsZero(N) then
			return MutableNullMat(dimP, M!.n, DefaultFieldOfMatrix(M!.gens[1]));
		elif not IsRModuleByRepresentationRep(N) then
			TryNextMethod();
		fi;

		phom := [ ];
		for i in Filtered([1..Size(lst)], k -> k in svec) do
			phom[i] := Hom(lst[i], M);
		od;
		hom0 := [ ];
		hom := [ ];
		for i in [1..Size(svec)] do
			px := Sum(List([1..i-1], q -> lst[svec[q]]!.n));
			for h in phom[svec[i]] do
				hx := MutableNullMat(dimP, M!.n, DefaultFieldOfMatrix(M!.gens[1]));
				hx{[px+1..px+Size(h)]}{[1..M!.n]} := h;
				Add(hom0, hx);
				hx_ := MutableNullMat(dimP, N!.n, DefaultFieldOfMatrix(N!.gens[1]));
				hx_{[px+1..px+Size(h)]}{[1..N!.n]} := h*psi;
				Add(hom, Flat(hx_));
			od;
		od;

		if IsZero(Size(hom)) then
			return fail;
		fi;
		if DefaultFieldOfMatrix(hom) = Rationals then
			coeff := SolutionIntMat(hom, Flat(phi));
		else
			coeff := SolutionMat(hom, Flat(phi));
		fi;
		if coeff = fail then
			return fail;
		fi;
		if IsRModuleOverZpOrderModp(N) then
			coeff := IntVecFFE(coeff);
		fi;
		ret := Sum([1..Size(hom)], i -> coeff[i]*hom0[i]);
		return ret;
	end
);

BlockIJ := function(A,i,j)
	local p, a;
	a := Filtered(A!.blockmap, x->x[1][1]=i and x[1][2]=j);
	if Size(a) = 0 then
		return [];
	else
		p := A!.blockidx[a[1][2]];
		return A!.gens{[p[1]..p[1]+p[2]-1]};
	fi;
end;

InstallMethod(ProjectiveIndecomposableForBasicOrder, [IsZpOrder and IsZpOrderMultiMatrixRep, IsPosInt],
	function(A, s)
		local B, M, Mrinv, rep, P, d;
		#if not (IsBound(A!.blockmap) and IsBound(A!.blockidx)) then
		#	TryNextMethod();
		#fi;
		if not IsBound(A!.gens_are_basis) then
			Print("# Warning! Generators not known to be a basis!\n");
		elif not A!.gens_are_basis then
			Print("# Warning! Generators not known to be a basis!\n");
		fi;

		if not s in [1..Size(SimpleModules(A))] then
			Error("Second argument not in range");
		fi;

		d := List(A!.idempot[s], v -> Filtered([1..Size(v)], q -> not IsZero(v[q][q])));
		return SubmoduleByBasisNC(DirectSumOfModules(Flat(List([1..Size(A!.nvec)],
			q -> List(d[q], z -> IrreducibleLattices(A)[q])))),
			Concatenation(List([1..Size(SimpleModules(A))],
				z -> List(BlockIJ(A, z, s), v -> Flat(List([1..Size(A!.nvec)], i -> v[i]{d[i]}))))));
	end
);

InstallMethod(ProjectiveIndecomposableLatticesForBasicOrder, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		if not IsBound(A!.pim) then
			A!.pim := List([1..Size(SimpleModules(A))], i -> ProjectiveIndecomposableForBasicOrder(A, i));
		fi;

		return List([1..Size(SimpleModules(A))], x -> [A!.pim[x], SimpleModules(A)[x]]);
	end
);

InstallMethod(ProjectiveCover, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M) # Returns [P, eta], where P is the projective cover of M and eta the epimorphism P -> M
		local A, pim, V, P0, epi, eta, phi, top;

		A := M!.order;
		if not IsBound(A!.pim) then
			Error(Concatenation("Projectives for the order of definition of M (first argument) have to be known!",
					" (Call ProjectiveIndecomposableLattices or ProjectiveIndecomposableLatticesForBasicOrder first!)"));
		fi;
		pim := A!.pim;
		V := TopEpimorphism(M);
		epi := List([1..Size(SimpleModules(A))], i -> HomToSimpleNC(pim[i], SimpleModules(A)[i])[1]);
		top := Flat(List([1..Size(SimpleModules(A))], i -> List([1..V[3][i]], j -> i)));
		P0 := DirectSumOfModules(List(top, i -> pim[i]));
		phi := DiagonalJoin(Concatenation(List([1..Size(SimpleModules(A))], i -> List([1..V[3][i]], j -> epi[i]))));
		eta := LiftHomomorphismToDirectSumNC(A!.pim, top, phi, M, V[2], V[1]);
		return [P0, eta];
	end
);

InstallMethod(ProjectiveCover, [IsRModuleOverZpOrder and IsZeroRModuleRep],
	function(M)
		return ZeroRModule(M!.order);
	end
);

InstallMethod(HellerTranslate, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M)
		local P, K;
		P := ProjectiveCover(M);
		if IsRModuleOverZpOrderModp(M) then
			K := BaseIntMat(Concatenation(List(NullspaceMat(P[2]), IntVecFFE), M!.p*IdentityMat(P[1]!.n)));
			return SubmoduleByBasisNC(P[1], K);
		elif IsRLatticeOverZpOrder(M) then
			return SubmoduleByBasisNC(P[1], NullspaceIntMat(P[2]));
		else
			TryNextMethod(); # Shouldn't happen...
		fi;
	end
);

InstallMethod(HellerTranslate, [IsRModuleOverZpOrder and IsZeroRModuleRep],
	function(M)
		return ZeroRModule(M!.order);
	end
);

InstallMethod(HellerTranslateModular, [IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep],
	function(M)
		local P, K;
		P := ProjectiveCover(M);
		K := NullspaceMat(P[2]);
		return SubmoduleByBasisNC(ReduceModP(P[1]), K);
	end
);

InstallMethod(HellerTranslateModular, [IsRModuleOverZpOrder and IsZeroRModuleRep],
	function(M)
		return ZeroRModule(M!.order);
	end
);

InstallMethod(GeneratorsForBasicOrder, [IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A)
		local B, i,j, blockmap, blockidx, rad, R, C, M, Mrinv, lst, rr, V, W, k, lincomb, gen, rad0, rad2,
		       rr1, rr2, t1, t2, lincomb_mat;
		if not (IsBound(A!.blockmap) and IsBound(A!.blockidx)) then
			TryNextMethod();
		fi;
		if not IsBound(A!.gens_are_basis) then
			Print("# Warning! Generators not known to be a basis!\n");
		elif not A!.gens_are_basis then
			Print("# Warning! Generators not known to be a basis!\n");
		fi;

		# Calculate Jac(e_i A e_i)
		rad := [ ];
		for i in [1..Size(SimpleModules(A))] do
			R := BlockIJ(A, i, i);
			# Regular representation:
			M := List(R, x -> Flat(x));
			Mrinv := TransposedMat(List(IdentityMat(Size(M)), v -> SolutionMat(TransposedMat(M), v)));
			lst := List(R, x -> [List(R, y -> Flat(List([1..Size(A!.nvec)],
				z -> y[z]*x[z]))) * Mrinv]);
			lst := Gcd(List(Flat(lst), DenominatorRat))*lst;
			C := ZpOrderByMultiMatrices(A!.p,lst);
			rr := RadicalOfModule(IrreducibleLattices(C)[1])[1]!.embedding_into_irr_lat[1];
			Add(rad, List(rr, x -> Sum(List([1..Size(x)], z -> x[z] * R[z]))));
		od;

		blockmap := [ ];
		blockidx := [ ];
		gen := [ ];
		for i in [1..Size(SimpleModules(A))] do
			for j in [1..Size(SimpleModules(A))] do
				lst := [];

				if IsZero(Size(BlockIJ(A, i, j))) then
					continue;
				fi;

				if i = j then
					# Calculate generators of e_i A e_i / Jac(e_i A e_i)
					R := BlockIJ(A,i,i);
					V := VectorSpace(GF(A!.p), List(rad[i], z -> One(GF(A!.p))*Flat(z)));
					for k in [1..Size(R)] do
						W := VectorSpace(GF(A!.p), [One(GF(A!.p))*Flat(R[k])]);
						if not IsSubspace(V, W) then
 							Add(lst, R[k]);
							V := V + W;
						fi;
					od;
					rad0 := rad[i];
				else
					rad0 := BlockIJ(A,i,j);
				fi;

				# Calculate e_i \Jac(A)^2 e_j	+ p*e_i*A*e_j
				rad2 := [ ];
				for k in [1..Size(SimpleModules(A))] do
					if k <> i then rr1 := BlockIJ(A,i,k); else rr1 := rad[k]; fi;
					if k <> j then rr2 := BlockIJ(A,k,j); else rr2 := rad[k]; fi;
					for t1 in rr1 do
						for t2 in rr2 do
							Add(rad2, List([1..Size(A!.nvec)], z -> t2[z]*t1[z]));
						od;
					od;
				od;
				rad2 := Concatenation(rad2, A!.p*BlockIJ(A, i, j));

				# Calculate a basis of e_i \Jac(A) e_j / e_i \Jac(A)^2 e_j (+ p*A)
				M := List(rad0, Flat);
				Mrinv := TransposedMat(List(IdentityMat(Size(M)), v -> SolutionMat(TransposedMat(M), v)));
				V := VectorSpace(GF(A!.p), List(rad2, z -> One(GF(A!.p))*(Flat(z)*Mrinv)));
				for k in [1..Size(M)] do
					W := VectorSpace(GF(A!.p), [One(GF(A!.p))*IdentityMat(Size(M))[k]]);
					if not IsSubspace(V, W) then
						Add(lst, rad0[k]);
						V := V + W;
					fi;
				od;

				if Size(lst) <> 0 then
					Add(blockidx, [Size(gen)+1, Size(lst)]);
					Add(blockmap, [[i, j], Size(blockidx)]);
					Append(gen, lst);
				fi;
			od;
		od;

		B := ZpOrderByMultiMatrices(A!.p, gen);
		if IsBound(A!.Qp_end_bases) then
			B!.Qp_end_bases := A!.Qp_end_bases;
		fi;
		if IsBound(A!.Qp_end_bases_smallmats) then
			B!.Qp_end_bases_smallmats := A!.Qp_end_bases_smallmats;
		fi;
		if IsBound(A!.component_names) then
			B!.component_names := A!.component_names;
		fi;
		if IsBound(A!.simple_names) then
			B!.simple_names := A!.simple_names;
		fi;
		if IsBound(A!.idempot) then
			B!.idempot := A!.idempot;
		fi;
		B!.blockmap := blockmap;
		B!.blockidx := blockidx;

		M := List(A!.gens, Flat);
		Mrinv := TransposedMat(List(IdentityMat(Size(M)), v -> SolutionMat(TransposedMat(M), v)));
		lincomb_mat := List(B!.gens, Flat) * Mrinv;
		lincomb := function(M)
			local N;
			N := RModuleOverZpOrder(B, List(lincomb_mat, z -> Sum(List([1..Size(z)], zz -> z[zz]*M!.gens[zz]))));
			if IsBound(M!.embedding_into_irr_lat) then
				N!.embedding_into_irr_lat := M!.embedding_into_irr_lat;
			fi;
			return N;
		end;

		B!.simple := List(A!.simple, z -> lincomb(z));
		B!.simple_multiplicities := A!.simple_multiplicities;

		return [B, lincomb];
	end
);

InstallGlobalFunction(HomForLattices, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep,
                                       IsRLatticeOverZpOrder and IsRModuleByRepresentationRep, IsList],
	function(M, N, opt)
		# opt: = [1], [2] or [1,2]: (opt = [ ] is possible, but useless...)
		#      1: Return Hom(M, N) embedded in Hom_Zp(M, N)
		#      2: Return Hom(M, N) embedded in Hom_Qp(M(x)Qp, N)x)Qp), using the matrices in
		#         Irr....!.Qp_end_basis_smallmats
		local irr, sourcedims, targetdims, hom0, i, j, k, l, k0, l0, m, n,
		      posx, posy, basa, alpha, beta, bas, bas_, smallbas, nf, det, trans, ret, lcm;

		if not (IsBound(M!.embedding_into_irr_lat) and IsBound(N!.embedding_into_irr_lat)) then
			Error("Embeddings into sums of irr. lattices must be known for source and target!");
		fi;
		if M!.order <> N!.order then
			Error("Modules defined over different orders!");
		fi;

		irr := IrreducibleLattices(M!.order);
		sourcedims := List(M!.embedding_into_irr_lat[2], k -> irr[k]!.n);
		targetdims := List(N!.embedding_into_irr_lat[2], k -> irr[k]!.n);
		alpha := M!.embedding_into_irr_lat[1];
		beta := _InverseRatMat@(N!.embedding_into_irr_lat[1]);
		det := Lcm(List(Flat(beta), DenominatorRat));
		beta := beta*det;

		# First construct the Qp basis of the Qp(x)Lambda-homomorphisms from Qp(x)M to Qp(x)N:
		bas := [ ];
		smallbas := [ ];
		for i in [1..Size(M!.embedding_into_irr_lat[2])] do
			for j in [1..Size(N!.embedding_into_irr_lat[2])] do
				posx := Sum(List([1..i-1], z -> irr[M!.embedding_into_irr_lat[2][z]]!.n));
				posy := Sum(List([1..j-1], z -> irr[N!.embedding_into_irr_lat[2][z]]!.n));
				if M!.embedding_into_irr_lat[2][i] = N!.embedding_into_irr_lat[2][j] then
					if not IsBound(irr[M!.embedding_into_irr_lat[2][i]]!.Qp_end_basis) then
						Error("Required endomorphism basis of irreducible lattice not known!");
					fi;
					if 2 in AsSet(opt) then
						if not IsBound(irr[M!.embedding_into_irr_lat[2][i]]!.Qp_end_basis_smallmats) then
							Error("Required endomorphism basis (in small rep.) of ireducible lattice not known!");
						fi;
						for n in irr[M!.embedding_into_irr_lat[2][i]]!.Qp_end_basis_smallmats do
							m := [ ];
							k := Size(Filtered(M!.embedding_into_irr_lat[2]{[1..i]},
							              z -> z = M!.embedding_into_irr_lat[2][i]));
							l := Size(Filtered(N!.embedding_into_irr_lat[2]{[1..j]},
							              z -> z = N!.embedding_into_irr_lat[2][j]));
							k0 := Size(Filtered(M!.embedding_into_irr_lat[2],
							              z -> z = M!.embedding_into_irr_lat[2][i]));
							l0 := Size(Filtered(N!.embedding_into_irr_lat[2],
							              z -> z = N!.embedding_into_irr_lat[2][j]));
							m[M!.embedding_into_irr_lat[2][i]] := MutableCopyMat(BlockMatrix([[k, l, n*det]],
							              k0, l0));
							Add(smallbas, m);
						od;
					fi;
					for n in irr[M!.embedding_into_irr_lat[2][i]]!.Qp_end_basis do
						m := List(alpha, v -> v{[1..Size(n)] + posx})*n*beta{[1..Size(n)] + posy};
						Append(bas, [Flat(m)]);
					od;
				fi;
			od;
		od;
		Unbind(alpha);
		Unbind(beta);

		# If bas is empty, there is only the zero homomorphism:
		if Size(bas) = 0 then
			if 1 in AsSet(opt) then Add(bas, [ ]); fi;
			if 2 in AsSet(opt) then Add(bas, [ ]); fi;
			return bas;
		fi;

		# Now intersect with Z:
		bas_ := Filtered(TransposedMat(bas), v -> not IsZero(v));
#		Print("Dim: ", Size(bas), " x ", Size(bas[1]), "   vs.      ", Size(bas_), " x ", Size(bas_[1]), "\n");
		bas_ := TransposedMatMutable(HermiteNormalFormIntegerMat(BaseIntMat(bas_)));
#		nf := NormalFormIntMat(bas_, 1 + 4 + 16);
#		lcm := Lcm(List([1..Size(nf.normal)], z -> nf.normal[z][z]));
#		trans := DiagonalMat(List([1..Size(nf.normal)], z -> nf.normal[z][z]^-1)) * nf.rowtrans;
#		trans := 1/lcm*HermiteNormalFormIntegerMat(lcm*trans);
		########### NEW ###############
	#	nf :=  SmithIntMatLLLTrans(bas_);
	#	lcm := Lcm(List([1..Size(nf[1])], z -> nf[1][z][z]));
	#	trans := DiagonalMat(List([1..Size(nf[1])], z -> nf[1][z][z]^-1)) * nf[2];
	#	trans := 1/lcm*HermiteNormalFormIntegerMat(lcm*trans);
		trans := _InverseRatMat@(bas_);
		########## END #################

		ret := [ ];
		if 1 in AsSet(opt) then
			bas := trans * bas;
			bas := List(bas, z -> List([1..Sum(sourcedims)],
		          zz -> z{[1..Sum(targetdims)] + (zz - 1)*Sum(targetdims)}));
			Add(ret, bas);
		fi;
		if 2 in AsSet(opt) then
			smallbas := List(trans, z -> Sum(List([1..Size(z)], zz -> z[zz]*smallbas[zz])));
			Add(ret, smallbas);
		fi;

		return ret;
	end
);

InstallMethod(HomToSimpleNC, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep,
                               IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep],
	function(M, N)
		local Mtr, Ntr, hom;

		Mtr := GModuleByMats(List(M!.gens, g -> TransposedMat(One(GF(M!.p))*g)), GF(M!.p));
		Ntr := GModuleByMats(List(N!.gens, g -> TransposedMat(g)), GF(N!.p));
		hom := MTX.Homomorphisms(Ntr, Mtr);

		return List(hom, m -> TransposedMatImmutable(m));
	end
);

InstallMethod(HomToSimpleNC, [IsRModuleOverZpOrder and IsZeroRModuleRep,
                               IsRModuleOverZpOrderModp and IsRModuleByRepresentationRep],
	function(M, N)
		return [ ];
	end
);

InstallMethod(Hom, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep,
                    IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(M, N)
		if M!.order <> N!.order then
			Error("Modules defined over different orders!");
		fi;
		if not (IsBound(M!.embedding_into_irr_lat) and IsBound(N!.embedding_into_irr_lat)
			and IsBound(M!.order!.Qp_end_bases)) then
			TryNextMethod();
		else
			return HomForLattices(M, N, [1])[1];
		fi;
	end
);

InstallMethod(Hom, [IsRModuleOverZpOrder and IsRModuleByRepresentationRep,
                    IsRModuleOverZpOrder and IsRModuleByRepresentationRep],
	function(M, N)
		local eqn, ns, nf, hom, i, j, k, s, t, x,
		      deltais, deltajt, range, hom0, R, modulus;
		# Check that M and N are defined over the same order:
		if M!.order <> N!.order then
			Error("Modules defined over different orders!");
		fi;

		# Check what we have to do (is N a torsion-module or a lattice?)
		if IsRModuleOverZpOrderModp(N) then
			modulus := N!.p;
			R := GF(modulus);
		elif IsRLatticeOverZpOrder(N) and IsRLatticeOverZpOrder(M) then
			modulus := 0;
			R := Integers;
		else # M torsion, N lattice => no homomorphisms
			return [ ];
		fi;

		# Construct the matrix eqn for the corresponding system of linear equations (over R):
		eqn := MutableNullMat(M!.n*N!.n, (M!.n*N!.n)*Size(M!.gens), R);
		for k in [1..Size(M!.gens)] do
			for i in [1..M!.n] do
				for j in [1..N!.n] do
					for s in [1..M!.n] do
						if i <> s then
							range := [j];
							deltais := 0;
						else
							range := [1..N!.n];
							deltais := 1;
						fi;
						for t in range do
							# Calculate (DeltaM(k)*e[i,j]-e[i,j]*DeltaN(k))[s,t]
							if j = t then
								deltajt := 1;
							else
								deltajt := 0;
							fi;
							eqn[i+M!.n*(j-1)][s + M!.n*(t-1) + M!.n*N!.n*(k-1)] :=
							   (M!.gens[k][s][i]*deltajt-N!.gens[k][j][t]*deltais) * One(R);
						od;
					od;
				od;
			od;
		od;

		# Calculate the kernel ns of eqn.
		if modulus = 0 then
			nf := NormalFormIntMat(eqn, 0+2+4+16);
			ns := nf.rowtrans;
			range := [nf.rank+1..Size(ns)];
		else # modulus = p
			ns := NullspaceMat(eqn);
			range := [1..Size(ns)];
		fi;

		# Reconstruct the homomorphisms from ns:
		hom := [];
		for x in range do
			hom0 := MutableNullMat(M!.n, N!.n, R);
			for i in [1..M!.n] do
				for j in [1..N!.n] do
					hom0[i][j] := ns[x][i+M!.n*(j-1)];
				od;
			od;
			Append(hom, [hom0]);
		od;

		return hom;
	end
);

InstallMethod(Hom, [IsRModuleOverZpOrder and IsZeroRModuleRep,
                    IsRModuleOverZpOrder],
	function(M, N)
		return [ ];
	end
);

InstallMethod(Hom, [IsRModuleOverZpOrder,
                    IsRModuleOverZpOrder and IsZeroRModuleRep],
	function(M, N)
		return [ ];
	end
);

InstallMethod(CondensationData, [IsGroup, IsGroup, IsCharacter],
	function(G, K, chi)
		local QK, emb, e, T, N, gens, epi;
		if not IsSubgroup(G, K) then
			Error("Second argument is supposed to be a subgroup of the first argument!");
		fi;
		if UnderlyingGroup(chi) <> K then
			Error("Third argument is supposed to be a character defined on the second argument");
		fi;
		QK := GroupRing(Rationals, K);
		emb := Embedding(K, QK);
		e := 1/Size(K)*Sum(List(K, k -> (k^-1)^chi*k^emb));
		N := Normalizer(G, K);
		T := InertiaSubgroup(N, chi);
		epi := NaturalHomomorphismByNormalSubgroup(T, K);
		gens := List(SmallGeneratingSet(Range(epi)), g -> PreImagesRepresentative(epi, g));
		Append(gens, Filtered(List(DoubleCosets(G, T, T), Representative), g -> not IsOne(g)));
		return rec( idempot := e, gens := gens, primes := AsSet(FactorsInt(Size(K))), G := G, K := K );
	end
);

InstallMethod(CondensationData, [IsGroup, IsGroup],
	function(G, K)
		return CondensationData(G, K, TrivialCharacter(K));
	end
);

InstallGlobalFunction(CondenseMatricesWithEvalMapNC, [IsFunction, IsRecord],
	function(hom, data)
		local emat, base, baseR, cgens;
		emat := CoefficientsAndMagmaElements(data.idempot);
		emat := Size(data.K)*Sum(List([1..Size(emat)/2], i -> hom(emat[2*i-1])*emat[2*i]));
		if IsZero(emat) then
			return [List([1..Size(data.gens)], i -> [[0]]), , ];
		fi;
		base := BaseIntMat(emat);
		baseR := RightInverse(base);
		cgens := List(data.gens, g -> base*emat*(hom(g))*emat*baseR);
		return [cgens, base, baseR];
	end
);

InstallGlobalFunction(CondenseMatricesNC, [IsList, IsList, IsRecord],
	function(Ggens, Greps, data)
		local hom;
		hom := GroupHomomorphismByImagesNC(data.G, GL(Size(Greps[1]), Integers), Ggens, Greps);
		return CondenseMatricesWithEvalMapNC(g -> g^hom, data);
	end
);

InstallGlobalFunction(CondenseTorsionRepNC, [IsList, IsList, IsRecord],
	function(Ggens, Greps, p, data)
		local hom, emat, base, baseR, cgens;
		hom := GroupHomomorphismByImagesNC(data.G, GL(Size(Greps[1]), GF(p)), Ggens, Greps);
		emat := CoefficientsAndMagmaElements(data.idempot);
		emat := Size(data.K)*Sum(List([1..Size(emat)/2], i -> emat[2*i-1]^hom*emat[2*i]));
		if IsZero(emat) then
			return List([1..Size(data.gens)], i -> [[Zero(GF(p))]]);
		fi;
		base := BaseMat(emat);
		baseR := RightInverse(base);
		cgens := List(data.gens, g -> base*emat*(g^hom)*emat*baseR);
		return cgens;
	end
);

InstallMethod(CondenseGroupRingNC, [IsZpOrder and IsZpOrderMultiMatrixRep,
                                      IsList, IsRecord],
	function(A, gens, data)
		local cdat, crep, nonzero, nonzero_simple, B, k, cgens;
		if IsBound(A!.eval_map) then
			cdat := List([1..Size(A!.nvec)], i -> CondenseMatricesWithEvalMapNC(A!.eval_map[i], data));
		else
			cdat := List([1..Size(A!.nvec)], i -> CondenseMatricesNC(gens, List(A!.gens, v -> v[i]), data));
		fi;
		crep := List(cdat, v -> v[1]);
		nonzero := Filtered([1..Size(crep)], v -> not ForAll(crep[v], vv -> IsZero(vv)));
		B := ZpOrderByMultiMatrices(A!.p, List([1..Size(data.gens)],
			i -> List(nonzero, n -> crep[n][i])));
		if IsBound(A!.Qp_end_bases) then
			B!.Qp_end_bases := List(nonzero, n -> List(A!.Qp_end_bases[n], e -> cdat[n][2]*e*cdat[n][3]));
		fi;
		if IsBound(A!.Qp_end_bases_smallmats) then
			B!.Qp_end_bases_smallmats := A!.Qp_end_bases_smallmats;
		fi;
		if IsBound(A!.component_names) then
			B!.component_names := A!.component_names{nonzero};
		fi;
		if IsBound(A!.simple) then # Condense the simple A-modules
			B!.simple := [ ];
			nonzero_simple := [ ];
			for k in [1..Size(A!.simple)] do
				cgens := CondenseTorsionRepNC(gens, A!.simple[k]!.gens, A!.p, data);
				if ForAny(cgens, z -> not IsZero(z)) then
					Add(B!.simple, RModuleOverZpOrder(B, cgens));
					Add(nonzero_simple, k);
				fi;
				B!.simple_multiplicities := A!.simple_multiplicities{nonzero}{nonzero_simple};
				if IsBound(A!.simple_names) then
					B!.simple_names := A!.simple_names{nonzero_simple};
				fi;
			od;
		fi;
		return B;
 	end
);

InstallMethod(CondenseGroupRingNC, [IsZpOrder, IsRecord],
	function(A, data) # Assumes that GeneratorsOfGroup(data.G) correspond to the generators of A
		return CondenseGroupRingNC(A, GeneratorsOfGroup(data.G), data);
	end
);

InstallMethod(CondensationProperties, [IsZpOrder, IsList, IsRecord],
	function(A, gens, data)
		local S, M, nonzero, dims, i, hom, emat;
		S := SimpleModules(A);
		M := DecompositionMatrix(A);
		nonzero := [ ];
		dims := [ ];
		for i in [1..Size(S)] do
			hom := GroupHomomorphismByImagesNC(data.G, GL(S[i]!.n, GF(A!.p)), gens, S[i]!.gens);
			emat := CoefficientsAndMagmaElements(data.idempot);
			emat := Sum(List([1..Size(emat)/2], j -> emat[2*j-1]^hom*emat[2*j]));
			if not IsZero(emat) then
				Add(nonzero, i);
				Add(dims, Rank(emat));
			fi;
		od;
		Print("Decomposition matrix:\n");
		M := Filtered(M{[1..Size(M)]}{nonzero}, v -> not IsZero(v));
		Display(M);
		Print("Dimensions of simple modules:\n");
		Display(dims);
		Print("Dimensions of irreducible lattices:\n");
		Display(TransposedMat(M*TransposedMat([dims]))[1]);
		Print("Lost ", Size(S) - Size(nonzero), " simple modules\n");
		Print("Lost ", Size(A!.nvec) - Size(M), " Wedderburn-components\n");
 	end
);

InstallMethod(CondensationProperties, [IsZpOrder, IsRecord],
	function(A, data)
		CondensationProperties(A, GeneratorsOfGroup(data.G), data);
	end
);

InstallMethod(DiagonalJoin, [IsList],
	function(L)
		local i, j, m, posr, posc, M;

		if Size(L) = 0 then
			Error("DiagonalJoin expects a list of at least one matrix");
		elif ForAny(L, l -> not IsMatrix(l)) then
			Error("DiagonalJoin expects its arguments to be matrices");
		fi;

		M := MutableNullMat(Sum(List(L, l -> Size(l))), Sum(List(L, l -> Size(TransposedMat(l)))),
			DefaultField(L[1][1][1]));
		posr := 0;
		posc := 0;
		for m in L do
			for i in [1..Size(m)] do
				for j in [1..Size(TransposedMat(m))] do
					M[posr + i][posc + j] := m[i][j];
				od;
			od;
			posr := posr + Size(m);
			posc := posc + Size(TransposedMat(m));
		od;

		return M;
	end
);

InstallMethod(DiagonalJoin, [IsMatrix, IsMatrix],
	function(M, N)
		return DiagonalJoin([M, N]);
	end
);

InstallMethod(RightInverse, [IsMatrix],
	function(M)
		local  F, V, i, Q, T, Mtr, lst;
		F := DefaultFieldOfMatrix(M);
		V := SubspaceNC(FullRowSpace(F, Size(M)), []);
		Q := MutableNullMat(Size(M[1]), Size(M), F);
		Mtr := TransposedMat(M);
		lst := [ ];
		i := 1;
		while (Dimension(V) < Size(M)) and (i <= Size(M[1])) do
			if Mtr[i] in V then
				i := i + 1;
				continue;
			fi;
			Q[i][Dimension(V) + 1] := 1;
			Add(lst, i);
			V := V + SubspaceNC(FullRowSpace(F, Size(M)), [ Mtr[i] ]);
		od;
		if Dimension(V) < Size(M) then
			return fail;
		fi;
		T := TransposedMat(Mtr{lst});
		if F = Rationals then
			return Q*_InverseRatMat@(T);
		else
			return Q*(T^-1);
		fi;
	end
);

# Specific setup for the symmetric groups

InstallGlobalFunction(ZpSnWedderburnComponentsNC, [IsPrime, IsList],
	function(p, part)
		local A, n, rho;
		n := Sum(part[1]);
		rho := List(part, p -> NaturalSpechtRepresentation(p));
		A := ZpOrderByMultiMatrices(p, List(GeneratorsOfGroup(SymmetricGroup(n)),
	     	g -> List(rho, p -> g^p)));
		A!.Qp_end_bases := List(A!.nvec, v -> [IdentityMat(v, Integers)]);
		A!.Qp_end_bases_smallmats := List(A!.nvec, v -> [IdentityMat(1, Integers)]);
		A!.component_names := List(part, p -> PartitionAsString(p));
		A!.eval_map := List(rho, x -> function(g)  return g^x; end);
		NameSimpleModulesByDims(A);
		return A;
	end
);

InstallGlobalFunction(PartitionAsString, [IsList],
	function(p)
		local ret, p0;
		p0 := Collected(p);
		Sort(p0, function(a,b) return a[1] > b[1]; end);
		p0 := List(p0, function(x) if x[2] = 1 then return String(x[1]);
	         else return Concatenation(String(x[1]),"^{", String(x[2]), "}"); fi; end);
		ret := Concatenation("(", p0[1], Concatenation(List(p0{[2..Size(p0)]}, x -> Concatenation(",", x))), ")");
		return ret;
	end
);

InstallGlobalFunction(ZpSn,
	function(p, n)
		return ZpSnWedderburnComponentsNC(p, Reversed(Partitions(n)));
	end
);

# Discriminants etc.:

InstallMethod(Valuation, [IsInt, IsInt],
	function(n, p)
		local n0, v;
		if not IsPrime(p) then
			Error("Second argument is supposed to be prime.");
		fi;
		n0 := n; v := -1;
		repeat
			n0 := n0/p; v := v + 1;
		until not IsInt(n0);
		return v;
	end
);

InstallMethod(Valuation, [IsRat, IsInt],
	function(r, p) # The (exponent-)valuation nu_p(r)
		 return Valuation(NumeratorRat(r), p) - Valuation(DenominatorRat(r), p);
	end
);

InstallMethod(GramMatrixOfTrace, [IsZpOrder and IsZpOrderMultiMatrixRep, IsList],
	function(A, u)
		local G, i, j;

		if not IsList(u) then
			Error("Second argument is supposed to be a list of rationals");
		elif not ForAll(u, uu ->  uu in Rationals) then
			Error("Second argument is supposed to be a list of rationals");
		elif not Size(u) = Size(A!.nvec) then
			Error("First argument has not the right number of components!");
		elif not IsBound(A!.gens_are_basis) then
			Print("# WARNING: The generators of the first argument are not known to be a basis!");
		elif not A!.gens_are_basis then
			Error("The generators of the first argument must be a basis of the order for this to work!");
		fi;

		G := MutableNullMat(Size(A!.gens), Size(A!.gens), Integers);
		for i in [1..Size(A!.gens)] do
			for j in [i..Size(A!.gens)] do
				G[i][j] := Sum(List([1..Size(A!.nvec)], k -> u[k]*Trace(A!.gens[i][k]*A!.gens[j][k])));
				G[j][i] := G[i][j];
			od;
		od;

		return G;
	end
);

InstallMethod(GramMatrixOfTrace, [IsZpOrder and IsZpOrderMatrixRep, IsRat],
	function(A, u)
		local G, i, j;
		if not IsBound(A!.gens_are_basis) then
			Print("# WARNING: The generators of the first argument are not known to be a basis!\n");
		elif not A!.gens_are_basis then
			Error("The generators of the first argument must be a basis of the order for this to work!");
		fi;

		G := MutableNullMat(Size(A!.gens), Size(A!.gens), Integers);
		for i in [1..Size(A!.gens)] do
			for j in [1..Size(A!.gens)] do
				G[i][j] := u*Trace(A!.gens[i]*A!.gens[j]);
				G[j][i] := G[i][j];
			od;
		od;

		return G;
	end
);

# Finds all suborders B of A which are minimal w.r.t. B \supseteq B^#
InstallMethod(SelfdualSuborders, [IsZpOrder and IsZpOrderMultiMatrixRep, IsList],
	function(A, u)
		local G, Asharp, Q, M, X, Xrinv, hom, S, N, T, i, j, h, hh, ns, cont, nb, c,
		      nb2, nb3, t, mm, gr, ll, b1, b2, MM, denom, lat, kk, qq, tt, ret;

		G := GramMatrixOfTrace(A, u);
		Asharp := ZpOrderByMultiMatrices(A!.p, List(Lcm(List(Flat(G^(-1)), DenominatorRat))*G^(-1),
			v -> Sum(List([1..Size(v)], j -> v[j]*A!.gens[j]))));
		b2 := Concatenation([[List(A!.nvec, z -> IdentityMat(z))], Asharp!.gens]);
		repeat
			b1 := b2;
			b2 := List(BaseIntMat(List(List(b1, z -> List(b1, zz -> List([1..Size(A!.nvec)], i -> z[i]*zz[i]))), Flat)), v -> UnFlattenMultiMatrixNC(v, A!.nvec));
		until b1 = b2;
		Asharp!.gens := b1;
		Asharp!.gens_are_basis := true;
		G := GramMatrixOfTrace(Asharp, u)^(-1);
		MM := List(G, v -> Sum(List([1..Size(v)], j -> v[j]*b2[j])));
		denom := Lcm(List(Flat(MM), DenominatorRat));
		MM := MM*denom;
		G := List(MM, m1 -> List(MM, m2 -> Sum(List([1..Size(A!.nvec)], i -> Trace(u[i]*m1[i]*m2[i])))));

		X := List(MM, Flat);
		Xrinv := TransposedMat(List(IdentityMat(Size(X)), v -> SolutionMat(TransposedMat(X), v)));
		Q := ZpOrderByMultiMatrices(A!.p, Concatenation(
			List(Asharp!.gens, g -> [List(MM, gg -> Flat(List([1..Size(A!.nvec)], i -> gg[i]*g[i])))*Xrinv]),
			List(Asharp!.gens, g -> [List(MM, gg -> Flat(List([1..Size(A!.nvec)], i -> g[i]*gg[i])))*Xrinv])
		));
		M := IrreducibleLattices(Q)[1];
		S :=	SimpleModules(Q);

		# Descend to a lattice with minimal discriminant (w.r.t. the condition that is has to contain its dual):
		N := M;
		T := IdentityMat(M!.n);
		while true do
			hom := List(S, s -> HomToSimpleNC(N, s));
			cont := false;
			for h in hom do
				for i in [0..Size(h)-1] do
					for j in Elements(FullRowSpace(GF(A!.p), i)) do
						hh := h[Size(h)-i] + Sum(List([1..i], z -> h[Size(h)-i+z]*j[z]));
						ns := NullspaceMat(hh);
						ns := BaseIntMat(Concatenation(List(ns, IntVecFFE), A!.p*IdentityMat(M!.n)));
						if Valuation(Lcm(List(Flat((ns*T*G*TransposedMat(T)*TransposedMat(ns))^(-1)), DenominatorRat)), A!.p) = 0 then
							T := ns*T;
							N := SubmoduleByBasisNC(N, ns);
							cont := true;
							break;
						fi;
					od;
					if cont then break; fi;
				od;
				if cont then break; fi;
			od;
			_DebugWriteByte@(STDOut, 43);
			if not cont then break; fi;
		od;
		_Debug@("  Found minimal lattice\n");

		# Calculate neighbors:
		nb := [HermiteNormalFormIntegerMat(T)];
		nb2 := [ ];
		while Size(nb) <> 0 do
			nb3 := [ ];
			for t in nb do
				gr := (t*G*TransposedMat(t))^(-1);
				mm := Filtered(List(MaximalSubmoduleBases(SubmoduleByBasisNC(M, t)), x -> TransposedMat(x[1])^(-1)*gr*t), z -> Valuation(Lcm(List(Flat(z), DenominatorRat)), A!.p) = 0);
				mm := List(mm, v -> v*Lcm(List(Flat(v), DenominatorRat))); # Make mm integral
				_Debug@(Size(mm), "\n");
				for j in [1..Size(mm)] do # Find unique representative of p-index
					c := NormalFormIntMat(mm[j], 1+4+8);
					mm[j] := c.rowtrans^(-1)*DiagonalMat(List([1..M!.n], z -> A!.p^Valuation(c.normal[z][z], A!.p)))*c.coltrans^(-1);
				od;
				for j in mm do
					N := SubmoduleByBasisNC(M, j);
					kk := MaximalSubmoduleBases(N);
					ll := [ ];
					for qq in kk do
						tt := (qq[1]*j)*G*TransposedMat(qq[1]*j);
						denom := Lcm(List(Flat(tt), DenominatorRat));
						if Maximum(List(ElementaryDivisorsMatDestructive(Integers, denom*tt),
							     z -> Valuation(DenominatorRat(1/z/denom), A!.p))) = 0 then
							Add(ll, qq);
						fi;
					od;
					Append(nb3, Filtered(List(ll, x -> HermiteNormalFormIntegerMat(x[1]*j)), z -> not z in nb and not z in nb2 and not z in nb3));
					_DebugWriteByte@(STDOut, 43);
				od;
				_Debug@(Size(nb3), " new lattices\n");
			od;
			Append(nb2, nb);
			nb := nb3;
		od;

		lat := List(nb2, nn -> List(nn, v -> Sum(List([1..Size(v)], i -> MM[i]*v[i]))));
		for j in [1..Size(lat)] do # Find unique representative of p-index
			c := NormalFormIntMat(List(lat[j], Flat), 1+4+8);
			lat[j] := List(
				HermiteNormalFormIntegerMat(c.rowtrans^(-1)*DiagonalMat(List([1..Size(c.normal)], z -> A!.p^Valuation(c.normal[z][z], A!.p)))*c.coltrans^(-1)), zz -> UnFlattenMultiMatrixNC(zz, A!.nvec));
		od;

		# Check which elements of lat are orders (uses the fact that 1 is contained in all lattices here):
		ret := [ ];
		for i in lat do
			kk := Product(ElementaryDivisorsMat(List(i, Flat)));
			qq := Product(ElementaryDivisorsMat(BaseIntMat(Concatenation(List(i, z ->
				List(i, zz ->Flat( List([1..Size(A!.nvec)], zzz -> z[zzz]*zz[zzz]))))))));
			if Valuation(kk, A!.p) = Valuation(qq, A!.p) then
				Add(ret, ZpOrderByMultiMatrices(A!.p, i));
			fi;
		od;

		return ret;
	end
);

InstallMethod(AreConjugate, [IsZpOrder and IsZpOrderMultiMatrixRep, IsZpOrder and IsZpOrderMultiMatrixRep],
	function(A, B)
		local Aop, C, V, i, j, L, j_, s, X, isgood, N, thislayer, lastlayer, old, max, m, mm, ll, irr;
		if A!.p <> B!.p then
			Error("Orders must be defined over the same ring");
		fi;

		if A!.nvec <> B!.nvec then
			return false;
		fi;

		Aop := ZpOrderByMultiMatrices(A!.p, List(A!.gens, x -> List(x, xx -> TransposedMat(xx))));
		C := ZpOrderByMultiMatrices(A!.p, Concatenation(
			List(Aop!.gens, x -> List(x, xx -> KroneckerProduct(xx, IdentityMat(Size(xx))))),
			List(B!.gens, x -> List(x, xx -> KroneckerProduct(IdentityMat(Size(xx)), xx)))
		));

		V := List([1..Size(C!.nvec)], z -> [ ]);
		for i in [1..Size(C!.nvec)] do
			irr := IrreducibleLattices(C);
			for j in [1..Size(irr)] do
				Unbind(irr[j]!.embedding_into_irr_lat);
			od;
			L := List(LatticeAlgorithm(irr[i], false), z -> SubmoduleByBasisNC(IrreducibleLattices(C)[i], z));
			for j in L do
				_DebugWriteByte@(STDOut, 45);
				isgood := true;
				j_ := RModuleOverZpOrder(Aop, j!.gens{[1..Size(Aop!.gens)]});
				for s in Filtered([1..Size(SimpleModules(Aop))], z -> DecompositionMatrix(Aop)[i][z] <> 0) do
					if Size(HomToSimpleNC(j_, SimpleModules(Aop)[s])) <> SimpleModules(Aop)[s]!.n then
						isgood := false;
						break;
					fi;
				od;
				if not isgood then
					continue;
				fi;
				j_ := RModuleOverZpOrder(B, j!.gens{Size(Aop!.gens) + [1..Size(B!.gens)]});
				for s in Filtered([1..Size(SimpleModules(B))], z -> DecompositionMatrix(B)[i][z] <> 0) do
					if Size(HomToSimpleNC(j_, SimpleModules(B)[s])) <> SimpleModules(B)[s]!.n then
						isgood := false;
						break;
					fi;
				od;
				if isgood then
					_Debug@("found one for ", i, "\n");
					if Size(irr) = 1 then
						return true;
					fi;
					Add(V[i], j);
				fi;
			od;
			_Debug@("done with ", i, "\n");
		od;

		for X in Cartesian(V) do
			_Debug@("====\n");
			N := DirectSumOfModules(X); # TODO: N has to be checked as well!
			lastlayer := [ IdentityMat(N!.n) ];
			thislayer := [ ]; old := [ ];
			while Size(lastlayer) <> 0 do
				_Debug@("New layer, looking at ", Size(lastlayer), " lattices\n");
				for m in lastlayer do
					max := List(MaximalSubmoduleBases(SubmoduleByBasisNC(N, m)), z -> HermiteNormalFormIntegerMat(z[1] * m));
					for mm in max do
						if mm in thislayer or mm in lastlayer or mm in old then
							continue;
						fi;
						ll := List([1..Size(X)], z -> List(mm, zz -> zz{Sum(List([1..z-1], kk -> X[kk]!.n)) + [1..X[z]!.n]}));
						ll := List(ll, z -> BaseMat(One(GF(A!.p)) * z));
						if ForAny([1..Size(X)], z -> Size(ll[z]) <> X[z]!.n) then
							continue;
						fi;
 						Add(thislayer, mm);

 						j := SubmoduleByBasisNC(N, mm);
 						isgood := true;
 						j_ := RModuleOverZpOrder(Aop, j!.gens{[1..Size(Aop!.gens)]});
 						for s in [1..Size(SimpleModules(Aop))] do
 							if Size(HomToSimpleNC(j_, SimpleModules(Aop)[s])) <> SimpleModules(Aop)[s]!.n then
 								isgood := false;
 								break;
 							fi;
 						od;
 						if not isgood then
 							continue;
 						fi;
 						j_ := RModuleOverZpOrder(B, j!.gens{Size(Aop!.gens) + [1..Size(B!.gens)]});
 						for s in [1..Size(SimpleModules(B))] do
 							if Size(HomToSimpleNC(j_, SimpleModules(B)[s])) <> SimpleModules(B)[s]!.n then
 								isgood := false;
 							break;
 						fi;
 						od;
 						if isgood then
 							return true;
 						fi;
					od;
				od;
				Append(old, lastlayer);
				lastlayer := thislayer;
				thislayer := [ ];
			od;
		od;

		return false;
	end
);

#################################### Newer functionality

# Returns integral matrices, but only a Q-basis of Hom_A(QM, QN). M and N need to be
# lattices over a group algebra, with generators being group elements (this is not checked).
InstallGlobalFunction(HomByReynoldsNC, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep,
                    IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(M, N)
    local G, Ggens, repM, repN, hom, rey, bas_, trans, reyG;
    Ggens := List([1..Size(M!.gens)], i -> DiagonalJoin(M!.gens[i], N!.gens[i]));
    G := GroupByGenerators(Ggens);
    repM := GroupGeneralMappingByImages(G, Ggens, M!.gens);
    repN := GroupGeneralMappingByImages(G, Ggens, N!.gens);
    rey := function(H)
				local max, M, f, cosets;
				if Size(H) = 1 then return X -> X; fi;
				max := MaximalSubgroupClassReps(H);
				M := Filtered(max, N -> Size(N) = Maximum(List(max, Size)))[1];
				cosets := List(RightCosetsNC(H, M), Representative);
				f := rey(M);
				return function(X)
					local X_;
					X_ := f(X);
					return Sum(List(cosets, g -> (g^(-1))^repM * X_ * g^repN));
				end;
    end;
    hom := List(Basis(MatrixSpace(Rationals, Dimension(M), Dimension(N))));
		reyG := rey(G);
    hom := List(hom, h -> reyG(h));
    hom := List(Basis(VectorSpace(Rationals, hom)));
    if Size(hom) > 0 then
        hom := List(Lcm(List(Flat(hom), DenominatorRat))*hom, Flat);
        bas_ := Filtered(TransposedMat(hom), v -> not IsZero(v));
        bas_ := TransposedMatMutable(HermiteNormalFormIntegerMat(BaseIntMat(bas_)));
        trans := _InverseRatMat@(bas_);
        hom := trans * hom;
        hom := List(hom, z -> List([1..Dimension(M)],  zz -> z{[1..Dimension(N)] + (zz - 1)*Dimension(N)}));
    fi;
    return hom;
	end
);

# Probabilistic test for isomorphism (used internally)
InstallGlobalFunction(_QuickIso@,
	function(L1, L2)
	  local hom12, i, h;
	  hom12 := HomForLattices(L1, L2, [1])[1];
	  for i in [1..30] do
		 h := Filtered(Sum(List(hom12, x -> Random(1, L1!.p)*x)), y -> true);
		 if Determinant(h) mod L1!.p <> 0 then
			return true;
		 fi;
	  od;
	  return false;
	end
);

# Some numerical invariants attached to a lattice (used internally)
InstallGlobalFunction(_LatticeInvariants@,
	function(L)
	  local endo, gram, inv, loewy, M, rad;

	  endo := List(HomForLattices(L, L, [2])[1], e -> Filtered(e, x -> true));
	  gram := List(endo, e1 -> List(endo, e2 -> Sum(List([1..Size(e1)], i -> Trace(e1[i]*e2[i])))));
	  inv := Filtered(ElementaryDivisorsMat(gram), x -> x <> 0);
	  inv := List(inv, x -> Valuation(x, L!.p));
	  Sort(inv);
	  M := ReduceModP(L);
	  loewy := [ ];
	  while Dimension(M) > 0 do
		 rad := RadicalOfModule(M);
		 M := rad[1];
		 Add(loewy, rad[2]);
	  od;
	  return [inv, loewy];
	end
);

InstallMethod(IsomorphismRModules, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep,
                    IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(L1, L2)
	  local end1, end2, hom12, hom21, c11, c22, p, A, B, V, bas, M, top, E, id1, id2, denom;

	  if L1!.n <> L2!.n then
		 return false;
	  fi;

	  p := L1!.p;
	  hom12 := List(HomForLattices(L1, L2, [2])[1], e -> Filtered(e, x ->true));
	  hom21 := List(HomForLattices(L2, L1, [2])[1], e -> Filtered(e, x ->true));

	  if Size(hom12) = 0 then # We assume that neither L1 nor L2 is the zero-module
		 return false;
	  fi;

	  end1 := Concatenation(List(hom12, alpha -> List(hom21, beta -> List([1..Size(alpha)],
		 i -> alpha[i]*beta[i]))));
	  denom := Lcm(List(Flat(end1), e -> DenominatorRat(e)));
	  id1 := List(end1[1], m -> IdentityMat(Size(m)));
	  if BaseIntersectionIntMats(denom*List(end1, Flat),
			[denom*Flat(id1)])[1][1] mod L1!.p = 0 then
		 return false;
	  fi;
	  end2 := Concatenation(List(hom21, alpha -> List(hom12, beta -> List([1..Size(alpha)],
		 i -> alpha[i]*beta[i]))));
	  denom := Lcm(List(Flat(end2), e -> DenominatorRat(e)));
	  id2 := List(end2[1], m -> IdentityMat(Size(m)));
	  if BaseIntersectionIntMats(denom*List(end2, Flat),
			[denom*Flat(id2)])[1][1] mod L1!.p = 0 then
		 return false;
	  fi;

	  hom12 := List(hom12, e -> DiagonalJoin(e));
	  hom21 := List(hom21, e -> DiagonalJoin(e));
	  end1 := List(end1, e -> DiagonalJoin(e));
	  end2 := List(end2, e -> DiagonalJoin(e));

	  # Check freeness of hom-spaces:
	  V := VectorSpace(Rationals, hom21);
	  bas := Basis(V, hom21);

	  # This is a somewhat dirty workaround, spinning would be nice
	  A := ZpOrderByMatrices(p, List(end1, e -> List(hom21, x -> Coefficients(bas, x*e))));

	  M := RModuleOverZpOrder(A, List(end1, e -> List(hom21, x -> Coefficients(bas, x*e))));
	  top := RadicalOfModule(M)[2];
	  E := List(SimpleModules(A), S -> Size(Hom(S, S)));
	  if not ForAll([1..Size(SimpleModules(A))], i -> SimpleModules(A)[i]!.n = top[i]*E[i]) then
		 return false;
	  fi;

	  V := VectorSpace(Rationals, hom12);
	  bas := Basis(V, hom12);
	  B := ZpOrderByMatrices(p, List(end2, e -> List(hom12, x -> Coefficients(bas, x*e))));
	  M := RModuleOverZpOrder(B, List(end2, e -> List(hom12, x -> Coefficients(bas, x*e))));
	  top := RadicalOfModule(M)[2];
	  E := List(SimpleModules(B), S -> Size(Hom(S, S)));
	  if not ForAll([1..Size(SimpleModules(B))], i -> SimpleModules(B)[i]!.n = top[i]*E[i]) then
		 return false;
	  fi;

	  return true;
	end
);

# Determines a subset of all sublattices of M, considering first M, then the
# maximal sublattices of M, then the maximal sublattices of the maximal sublattices,
# and so on. A lattice is dropped if it does not filfill the condition cond or if
# it is isomorphic to a lattice that found before. If a lattice is dismissed,
# we donot consider its sublattices. (this function is mostly for internal use)
InstallGlobalFunction(AllLatticesCond, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep, IsFunction],
	function(M, cond)
		local simple, lat, lastlayer, thislayer, m, n, x, max, bas, old, inv, cand, perm;
		simple := SimpleModules(M!.order);
		lat := [ ];

		lastlayer := [ [ M, _LatticeInvariants@(M) ] ];

		while Size(lastlayer) > 0 do
			_Debug@("New layer, found so far: ", Size(lat), "; looking at ", Size(lastlayer), " lattices", "\n");
			thislayer := [ ];

			for m in lastlayer do
				max := List(MaximalSubmodules(m[1]), x -> x[1]);

				Append(lat, [m]);
				# Test which lattices are isomorphic to ones we already found
				for n in max do
					# Print(Position(max, n), " of ", Size(max), "\n");
					_DebugWriteByte@(STDOut, 46);
					if not cond(n) then
					  continue;
					fi;

					inv := _LatticeInvariants@(n);
					cand := Filtered(Concatenation(thislayer, lat), x -> x[2] = inv);

					old := false;
					for x in cand do
					  if _QuickIso@(n, x[1]) then
					     old := true;
					     break;
					  fi;
					od;
					if old then
					  continue;
					fi;
					if ForAny(cand, x -> IsomorphismRModules(x[1], n)) then
						continue;
					fi;
					_DebugWriteByte@(STDOut, 33);
					Add(thislayer, [n, inv]);
				od;
			od;
			_Debug@("\n");
			lastlayer := thislayer;
		od;

		return List(lat, L -> L[1]);
	end
);

# Returns a list of all isoclasses of sublattices of M
InstallGlobalFunction(AllLattices, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(M)
  	return AllLatticesCond(M, x -> true);
	end
);

# Returns a list of all isoclasses sublattices of DirectSum(L1, L2) which surject onto
# both L1 and L2
InstallGlobalFunction(AllLatticesInd, [IsRLatticeOverZpOrder and IsRModuleByRepresentationRep, IsRLatticeOverZpOrder and IsRModuleByRepresentationRep],
	function(L1, L2)
	  local L, d1, d2, i1, i2, surj;
	  d1 := L1!.n; d2 := L2!.n;
	  i1 := Valuation(Determinant(L1!.embedding_into_irr_lat[1]), L1!.p);
	  i2 := Valuation(Determinant(L2!.embedding_into_irr_lat[1]), L2!.p);
	  L := DirectSumOfModules(L1, L2);
	  surj := function(M)
		 local v1, v2;
		 v1 := Determinant(BaseIntMat(List(M!.embedding_into_irr_lat[1], x -> x{[1..d1]})));
		 v2 := Determinant(BaseIntMat(List(M!.embedding_into_irr_lat[1], x -> x{[d1+1..d1+d2]})));
		 if Valuation(v1, L!.p) <> i1 or Valuation(v2, L!.p) <> i2 then
			return false;
		 else
			return true;
		 fi;
	  end;
	  return AllLatticesCond(L, surj);
	end
);

# This function takes a list "lst" of rational matrices which generate a
# Z-subalgebra of GL(n,Q) which is finitely generated as a Z-modue (not checked),
# and returns a list "ret" of matrices in GL(n, Z) such that ret[i] = T*lst[i]*T^(-1)
# for all i and some T in GL(n,Q).
# This function is included for convenience, as defining an order requires the representation
# to be integral. Input is not checked.
InstallGlobalFunction(SpinningAlgorithmNC,
	function(lst)
		local bas, d, ret;
		bas := IdentityMat(Size(lst[1]));
		repeat
			bas := Concatenation(bas,
				Concatenation(List(lst, g -> bas * g)));
			d := Lcm(List(Flat(bas), DenominatorRat));
			bas := 1/d * BaseIntMat(d*bas);
			ret := List(lst, g -> bas*g*bas^(-1));
		until ForAll(Flat(ret), IsInt);
		return ret;
	end
);

# Mostly for internal use, but may be useful elsewhere
InstallGlobalFunction(CoxeterLength, [IsPerm, IsPosInt],
	function(sigma, n)
		return Sum(List([1..n], i -> Size(Filtered([i+1..n], j -> i^sigma > j^sigma))));
	end
);

# Mostly for internal use, but may be useful elsewhere
InstallGlobalFunction(EnumerateStandardTableaux, [IsList],
	function(lambda)
		local recursive, lambda2;
		lambda2 := StructuralCopy(lambda);
		while Size(lambda2) <= Sum(lambda) do Add(lambda2, 0); od;

		recursive := function(lambda)
			local n, i, lambda_, lst, t, ret;
			n := Sum(lambda);
			if n = 1 then
				return [ [ [ 1 ], [ ] ] ];
			else
				ret := [ ];
				for i in [1..n] do
					if lambda[i] - 1 >= lambda[i + 1]  then
						lambda_ := StructuralCopy(lambda);
						lambda_[i] := lambda[i] - 1;
						lst := List(EnumerateStandardTableaux(lambda_), t -> Concatenation(t, [[ ]]));
						for t in lst do
							Add(t[i], n);
						od;
						Append(ret, lst);
					fi;
				od;
				return ret;
			fi;
		end;

		return List(recursive(lambda2), x -> Filtered(x, y -> Size(y) <> 0));
	end
);

# Mostly for internal use, but may be useful elsewhere
InstallGlobalFunction(ContentVector, [IsList],
	function(T)
		local a, i, j;
		a := List([1..Size(Flat(T))], i -> 0);
		for i in [1..Size(T)] do
			for j in [1..Size(T[i])] do
				a[T[i][j]] := j - i;
			od;
		od;
		return a;
	end
);

# Mostly for internal use, but may be useful elsewhere
InstallGlobalFunction(IsStdTableau, [IsList],
	function(T)
		if ForAll([2..Size(T)], i -> ForAll([1..Size(T[i])], j -> T[i][j] > T[i-1][j])) then
			if ForAll([1..Size(T)], i -> ForAll([2..Size(T[i])], j -> T[i][j] > T[i][j-1])) then
				return true;
			fi;
		fi;
		return false;
	end
);

# Implements Young's "natural" Specht representation for S^lambda, where lambda
# is a partition of some n >= 1. See G.D.James: "The Representation Theory of
# the Symmetric Group", Example 25.2
InstallGlobalFunction(NaturalSpechtRepresentation, [ IsList ],
	function(lambda)
		local n, j, i, i_, std, rho, ai, lst, thislayer, newlayer, thislayer_pos,
			newlayer_pos, Tj, pos, e, eI, gen;
		n := Sum(lambda);
		if n = 1 then
			return GroupHomomorphismByImages(SymmetricGroup(1), GL(1, Integers), [ ], [ ]);
		fi;
		std := EnumerateStandardTableaux(lambda);
		lst := [ ];
		for j in [1..n-1] do
			rho := MutableNullMat(Size(std), Size(std));
			for i in [1..Size(std)] do
				ai := ContentVector(std[i]);
				if ai[j + 1] = ai[j] + 1 then
					rho[i][i] := 1;
				elif ai[j + 1] = ai[j] - 1 then
					rho[i][i] := -1;
				else
					rho[i][i] := 1/(ai[j+1]-ai[j]);
					i_ := Position(std, List(std[i], x -> List(x, y -> y^(j,j+1))));
					if CoxeterLength(PermList(Flat(std[i])), n)
								< CoxeterLength(PermList(Flat(std[i_])), n) then
						rho[i][i_] := 1;
					else
						rho[i][i_] := 1 - 1/(ai[j+1]-ai[j])^2;
					fi;
				fi;
			od;
			Add(lst, rho);
		od;

		e := List([1..Size(std)], x -> 0);
		e[1] := IdentityMat(Size(std))[1];
		thislayer := [ std[1] ];
		thislayer_pos := [ 1 ];
		while Size(thislayer) > 0 do
			newlayer := [ ];
			newlayer_pos := [ ];
			for i in [1..Size(thislayer)] do
				for j in [1..n-1] do
					Tj := List(thislayer[i], x -> List(x, y -> y^(j, j+1)));
					if IsStdTableau(Tj) then
						pos := Position(std, Tj);
						if e[pos] = 0 then
							Add(newlayer, Tj);
							Add(newlayer_pos, pos);
							e[pos] := e[thislayer_pos[i]]*lst[j];
						fi;
					fi;
				od;
			od;
			thislayer := newlayer;
			thislayer_pos := newlayer_pos;
		od;
		eI := Inverse(e);
		lst := List(lst, x -> e*x*eI);
		gen := List([1..n-1], j -> (j, j+1));
		return GroupHomomorphismByImagesNC(SymmetricGroup(n), GL(Size(std), Integers),
				gen, lst);
	end
);

MakeReadWriteGlobal("_DEBUG_OUTPUT@");
_DEBUG_OUTPUT@ := false;

InstallGlobalFunction(SetDebugOutput, [IsBool],
	function(b)
		_DEBUG_OUTPUT@ := b;
	end
);

InstallGlobalFunction(_Debug@,
	function(arg)
		if _DEBUG_OUTPUT@ then
			CallFuncList(Print, arg);
		fi;
	end
);

InstallGlobalFunction(_DebugWriteByte@,
	function(arg)
		if _DEBUG_OUTPUT@ then
			CallFuncList(WriteByte, arg);
		fi;
	end
);

##### Pretty output

InstallGlobalFunction(CreateHTMLSummary, [IsString, IsString, IsString, IsZpOrder and IsZpOrderMultiMatrixRep],
	function(fd, title, gapfileprefix, A, loewyN)
		local PartitionLatex2Html, out, i, j, k, P, loewyMOD, loewyINT, xx;

		PartitionLatex2Html := function(s)
			return ReplacedString(ReplacedString(s, "^{", "<sup>"), "}", "</sup>");
		end;

		out := OutputTextFile(fd, false);
		SetPrintFormattingStatus(out, false);
		PrintTo(out, "<center><h1>", title, "</h1>\n");

		PrintTo(out, "<div><h2>Generators</h2>");
		PrintTo(out, "<i>Remark: For convenience there are two files: one file contains a generating set for the order, the other one an entire basis. The orders they generate are the same.</i><br/>");
		PrintTo(out, "<br/><a href=\"", gapfileprefix, "_gens.g\">Download GAP file (generators)</a><br/>");
		PrintTo(out, "<a href=\"", gapfileprefix, "_basis.g\">Download GAP file (basis)</a></div>");

		PrintTo(out, "<div><h2>Decomposition matrix</h2>");
		PrintTo(out, "<table border=\"0\" cellspacing=\"20\">\n");
		PrintTo(out, "<tr><td></td><td></td>");
		for i in [1..Size(SimpleModules(A))] do
			PrintTo(out, "<td><b>S<sub>", i, "</sub></b></td>");
		od;
		PrintTo(out, "</tr>");
		if IsBound(A!.simple_names) then
			PrintTo(out, "<tr><td></td><td></td>");
			for i in [1..Size(SimpleModules(A))] do
				PrintTo(out, "<td><b>", PartitionLatex2Html(A!.simple_names[i]), "</b></td>");
			od;
			PrintTo(out, "</tr>");
		fi;
		for j in [1..Size(A!.nvec)] do
			PrintTo(out, "<tr><td><b>A<sub>", j, "</sub> = M<sub>", A!.nvec[j], "</sub>(Z<sub>", A!.p, "</sub>)</b></td>");
			if IsBound(A!.component_names) then
				PrintTo(out, "<td><b>", PartitionLatex2Html(A!.component_names[j]), "</b></td>");
			else
				PrintTo(out, "<td></td>");
			fi;
			for i in [1..Size(SimpleModules(A))] do
				PrintTo(out, "<td>");
				if A!.simple_multiplicities[j][i] <> 0 then
					PrintTo(out, A!.simple_multiplicities[j][i]);
				else
					PrintTo(out, ".");
				fi;
				PrintTo(out, "</td>");
			od;
			PrintTo(out, "</tr>");
		od;
		PrintTo(out, "</table>");
		PrintTo(out, "</div>");

		xx := GeneratorsForBasicOrder(A);
		P := List([1..Size(SimpleModules(A))], i -> xx[2](ProjectiveIndecomposableForBasicOrder(A, i)));
		loewyINT := List(P, p -> RadicalSeries(p, loewyN));
		loewyMOD := List(P, p -> RadicalSeries(ReduceModP(p), loewyN));
		PrintTo(out, "<div><h2>Loewy layers</h2>");
		for k in [1..Size(SimpleModules(A))] do
			PrintTo(out, "<h3>Projective cover of S<sub>", k, "</sub></h3>\n");
			PrintTo(out, "<table border=\"0\" cellspacing=\"20\">\n<tr><td></td><td><b>Modular part</b></td><td><b>Integral part</b></td></tr>");
			for i in [1..Size(loewyINT[k])] do
				PrintTo(out, "<tr><td><b>Layer ", i, "</b></td><td>");
				for j in [1..Size(SimpleModules(A))] do
					if loewyMOD[k][i][j] > 0 then
						PrintTo(out, "S<sub>", j, "</sub>");
						if loewyMOD[k][i][j] > 1 then
							PrintTo(out, "<sup>", loewyMOD[k][i][j], "</sup>");
						fi;
						if ForAny(loewyMOD[k][i]{[j+1..Size(SimpleModules(A))]}, v -> v > 0) then
							PrintTo(out, ", ");
						fi;
					fi;
				od;
				PrintTo(out, "</td><td>");
				for j in [1..Size(SimpleModules(A))] do
					if loewyINT[k][i][j] - loewyMOD[k][i][j] > 0 then
						PrintTo(out, "S<sub>", j, "</sub>");
						if loewyINT[k][i][j] - loewyMOD[k][i][j] > 1 then
							PrintTo(out, "<sup>", loewyINT[k][i][j] - loewyMOD[k][i][j], "</sup>");
						fi;
						if ForAny([j+1..Size(SimpleModules(A))], v ->  loewyINT[k][i][v] - loewyMOD[k][i][v] > 0) then
							PrintTo(out, ", ");
						fi;
					fi;
				od;
				PrintTo(out, "</td>");
			od;
			PrintTo(out, "</table>");
		od;
		PrintTo(out, "</div>");

		PrintTo(out, "</center>");
		CloseStream(out);
	end
);

InstallGlobalFunction(PrintDecompositionMatrixAsLatex,
	function(A, stream...)
		local i, j, p, s0, t0, out, dec;

		if Length(stream) = 0 then
			out := STDOut;
		elif Length(stream) = 1 then
			out := stream[1];
		else
			Error("Too many arguments!");
		fi;

		dec := DecompositionMatrix(A); # We need to call this to make sure the revelant data is present
		PrintTo(out, "$$");
		if IsBound(A!.component_names) then
			s0 := "cc"; t0 := "&";
		else
			s0 := "c"; t0 := "";
		fi;
		PrintTo(out, "\\begin{array}{||", s0, "|", Concatenation(List([1..Size(A!.simple)], i -> "c")), "||} \\hline\n");
		if IsBound(A!.simple_names) then
			PrintTo(out, t0, Concatenation(List([1..Size(A!.simple)], i -> Concatenation(" & ", A!.simple_names[i]))), " \\\\");
		fi;
		PrintTo(out, t0, Concatenation(List([1..Size(A!.simple)], i -> Concatenation(String(" & "),
		 	String(Dimension(SimpleModules(A)[i]))))), " \\\\ \\hline");
		for i in [1..Size(A!.nvec)] do
			if IsBound(A!.component_names) then
				PrintTo(out, A!.component_names[i], "&\n");
			fi;
			PrintTo(out, A!.nvec[i]);
			for j in [1..Size(A!.simple)] do
				p := A!.simple_multiplicities[i][j];
				if IsZero(p) then p := "."; fi;
				PrintTo(out, " & ", p, "\n");
			od;
			PrintTo(out, "\\\\\n");
		od;
		PrintTo(out, "\\hline\\end{array}\n");
		PrintTo(out, "$$\n");
	end
);

InstallGlobalFunction(PrintBasisOfOrderAsMarkdown,
	function(A, stream...)
		local i, j, k, l, i0, j0,  pos, idx, supp, rowidx, colidx, p, s0, t0, out;

		if Length(stream) = 0 then
			out := STDOut;
		elif Length(stream) = 1 then
			out := stream[1];
		else
			Error("Too many arguments!");
		fi;

		for i in [1..Size(A!.simple)] do
			for j in [1..Size(A!.simple)] do
				pos := Filtered(A!.blockmap, x -> x[1] = [i, j]);
				if Size(pos) <> 0 then
	                PrintTo(out, "* $\\textrm{Hom}_\\Lambda(\\mathcal{P}(S_", i, "), \\mathcal{P}(S_", j , "))^{\\top}$:\n\n");
					idx := A!.blockidx[pos[1][2]];
					# supp := Filtered([1..Size(A!.nvec)], z -> not IsZero(A!.idempot[i][z]) and not IsZero(A!.idempot[j][z]));
					supp := [1..Size(A!.nvec)];
					colidx := List(supp, z -> Filtered([1..A!.nvec[z]], zz -> A!.idempot[i][z][zz][zz] <> 0));
					rowidx := List(supp, z -> Filtered([1..A!.nvec[z]], zz -> A!.idempot[j][z][zz][zz] <> 0));
					PrintTo(out, "$$\n \\left[");
					for l in [1..Size(supp)] do
	                    PrintTo(out, "\\left(\\begin{array}{", Concatenation(List([1..A!.nvec[supp[l]]], z -> "c")), "}\n");
						for i0 in [1..A!.nvec[supp[l]]] do
							for j0 in [1..A!.nvec[supp[l]]] do
								if i0 in rowidx[l] and j0 in colidx[l] then
									PrintTo(out, " \\scriptscriptstyle \\blacksquare \n");
								else
									PrintTo(out, " \\scriptscriptstyle \\square \n");
								fi;
								if j0 <> A!.nvec[supp[l]] then
									PrintTo(out, " & ");
								elif i0 <> A!.nvec[supp[l]] then
									PrintTo(out, " \\\\\n");
								fi;
							od;
						od;
	                    PrintTo(out, "\\end{array}\\right)");
						if l <> Size(supp) then
							PrintTo(out, ",\n");
						fi;
					od;
	                PrintTo(out, "\n \\right] $$\n<br/>");
					supp := Filtered([1..Size(A!.nvec)], z -> not IsZero(A!.idempot[i][z]) and not IsZero(A!.idempot[j][z])); # ...
					colidx := List(supp, z -> Filtered([1..A!.nvec[z]], zz -> A!.idempot[i][z][zz][zz] <> 0));
					rowidx := List(supp, z -> Filtered([1..A!.nvec[z]], zz -> A!.idempot[j][z][zz][zz] <> 0));
	                PrintTo(out, "$\\begin{array}{||", Concatenation(List(supp, z -> "c|")), "|}\n\\hline");
					for k in [idx[1]..idx[1]+idx[2]-1] do
						for l in [1..Size(supp)] do
							PrintTo(out, "\\begin{array}{", Concatenation(List(colidx[l], z -> "c")), "}\n");
							for i0 in [1..Size(rowidx[l])] do
								for j0 in [1..Size(colidx[l])] do
									p := A!.gens[k][supp[l]][rowidx[l][i0]][colidx[l][j0]];
									if IsZero(p) then p := "\\cdot"; fi;
									PrintTo(out, " ", p);
									if j0 <> Size(colidx[l]) then
										PrintTo(out, " &\n");
									elif i0 <> Size(rowidx[l]) then
										PrintTo(out, " \\\\\n");
									fi;
								od;
							od;
							PrintTo(out, "\n\\end{array}");
							if l <> Size(supp) then
								PrintTo(out, "\n&\n");
							fi;
						od;
						PrintTo(out, " \\\\\\hline\n");
					od;
					PrintTo(out, "\\end{array}$\n");
				fi;
			od;
		od;
	end
);

if IsBound(JupyterRender) then
	InstallGlobalFunction(JupyterDisplayDecompositionMatrix, [IsZpOrder and IsZpOrderMultiMatrixRep],
		function(A)
			local str, out;
			str := "";
			out := OutputTextString(str,true);
			PrintDecompositionMatrixAsLatex(A, out);
			return JupyterRenderable(rec(("text/latex") := str), rec());
		end
	);


	InstallGlobalFunction(JupyterDisplayBasisOfOrder, [IsZpOrder and IsZpOrderMultiMatrixRep],
		function(A)
			local str, out;
			str := "";
			out := OutputTextString(str,true);
			PrintBasisOfOrderAsMarkdown(A, out);
			return JupyterRenderable(rec(("text/markdown") := str), rec());
		end
	);

	InstallGlobalFunction(JupyterDisplayMultiMatrix, [IsList],
		function(lst)
			local str, n;
			str := "$\\left[";
			for n in [1..Size(lst)] do
				Append(str, LaTeX(lst[n]));
				if n <> Size(lst) then Append(str, ", "); fi;
			od;
			Append(str, "\\right]$");
			return JupyterRenderable(rec(("text/latex") := str), rec());
		end
	);
fi;
