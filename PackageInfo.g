SetPackageInfo( rec(

PackageName := "orders",
Subtitle := "Orders over the p-adic integers and their modules",
Version := "1.0",
Date := "09-07-2018",

Persons := [

 rec(
      LastName      := "Eisele",
      FirstNames    := "Florian",
      IsAuthor      := true,
      IsMaintainer  := true,
      Email         := "florian.eisele@city.ac.uk",
      WWWHome       := "",
      PostalAddress := Concatenation("City, University of London\n",
	                                "Northampton Square\n",
	                                "London EC1V 0HB\n",
	                                "United Kingdom"),
      Place         := "London",
      Institution   := "City, University of London"),

],

Status := "none",

README_URL := "",
PackageInfoURL := "",

AbstractHTML :=
"",
PackageWWWHome := "",

PackageDoc := rec(
  BookName  := "orders",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.htm",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Computing with orders over the p-adic integers",
  Autoload  := false),

Dependencies := rec(
  GAP := ">=4.8",
  NeededOtherPackages := [ ],
  SuggestedOtherPackages := [],
  ExternalConditions := [] ),

AvailabilityTest := ReturnTrue,
BannerString := Concatenation("--------------------------------------------------------------\n",
                              "  orders - Computing with orders over the p-adic integers  \n",
                              "                   (version ", ~.Version, ")\n\n",
                              "      Florian Eisele (florian.eisele@city.ac.uk)          \n",
                              "--------------------------------------------------------------\n"),
Autoload := false,
Keywords := ["orders", "lattices", "group rings"]

));
