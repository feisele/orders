SetPackageInfo( rec(

PackageName := "orders",
Subtitle := "Orders over the p-adic integers and their modules",
Version := "1.0",
Date := "14/12/2020",

Persons := [

 rec(
      LastName      := "Eisele",
      FirstNames    := "Florian",
      IsAuthor      := true,
      IsMaintainer  := true,
      Email         := "florian.eisele@city.ac.uk",
      WWWHome       := "https://feisele.github.io",
      PostalAddress := Concatenation("City, University of London\n",
	                                "Northampton Square\n",
	                                "London EC1V 0HB\n",
	                                "United Kingdom"),
      Place         := "London",
      Institution   := "City, University of London"),

],

Status := "other",

README_URL := "https://github.com/feisele/orders/blob/master/README.md",
PackageInfoURL := "https://github.com/feisele/orders/blob/master/README.md",

AbstractHTML :=
"",
PackageWWWHome := "https://github.com/feisele/orders/",

ArchiveURL := "https://github.com/feisele/orders/archive/1.0.tar.gz",

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
Keywords := ["orders", "lattices", "group rings"],
ArchiveFormats := ""

));
