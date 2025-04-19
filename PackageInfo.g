SetPackageInfo( rec(

PackageName := "orders",
Subtitle := "Orders over the p-adic integers and their modules",
Version := "1.0.3",
Date := "19/04/2025",
License := "GPL-2.0-or-later",

Persons := [

 rec(
      LastName      := "Eisele",
      FirstNames    := "Florian",
      IsAuthor      := true,
      IsMaintainer  := true,
      Email         := "florian.eisele@manchester.ac.uk",
      WWWHome       := "https://feisele.github.io",
      PostalAddress := "University of Manchester",
      Place         := "Manchester",
      Institution   := "University of Manchester"),

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
                              "      Florian Eisele (florian.eisele@manchester.ac.uk)          \n",
                              "--------------------------------------------------------------\n"),
Autoload := false,
Keywords := ["orders", "lattices", "group rings"],
ArchiveFormats := ""

));
