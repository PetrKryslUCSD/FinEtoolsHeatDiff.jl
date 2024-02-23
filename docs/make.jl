using Documenter, FinEtools, FinEtoolsHeatDiff

makedocs(
	modules = [FinEtoolsHeatDiff],
	warnonly = Documenter.except(:linkcheck, :footnote),
	doctest = false, clean = true,
	format = Documenter.HTML(), # prettyurls = false
	authors = "Petr Krysl",
	sitename = "FinEtoolsHeatDiff.jl",
	pages = Any[
	"Home" => "index.md",
	"Guide" => "guide/guide.md",
	"Types and Functions" => Any[
		"man/man.md"]
		]
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsHeatDiff.jl.git",
)
