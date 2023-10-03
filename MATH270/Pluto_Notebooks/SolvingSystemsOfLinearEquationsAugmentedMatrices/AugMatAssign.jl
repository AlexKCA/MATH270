### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 2fa9b038-5e2c-11ee-2a7f-29e9335e0135
using LinearAlgebra, Latexify, LaTeXStrings

# ╔═╡ 076ea2b4-fcc0-4dbe-8242-ef7344e982f0
function sysEq(n)
	# Returns a system of n linear equations in n unknowns
	A = []
	strings = String[]
	detA = 0
	
	while  detA == 0
		A = rand(-9:9, (n,n))
		detA = det(A)
	end
	
	b = rand(-9:9, (n,1))
	Ab = [A b]
	
	for i in 1:n
		if  Ab[i,1] != 0 
				str = "$(Ab[i,1])x_1"
			else str = ""
		end
		for j in 2:n
			if Ab[i,j] < 0
				str *= "-$(-Ab[i,j])x_{$j}"
			elseif Ab[i,j] > 0
				str *= "+$(Ab[i,j])x_{$j}"
			end
		end
		str *= "&=$(Ab[i,n+1])"
		push!(strings,str)
	end
	str = ""
	for i in 1:(n)
		str *= "$(strings[i])\\\\"
	end
	eqnalign = "\\begin{align}" * str * "\\end{align}"
	return convert(Array{Rational{Int64},2},Ab), eqnalign
end

# ╔═╡ 3e5e4708-7e41-4d62-8f50-806f8d47fe67
function equations(n)
	return LaTeXString(sysEq(n)[2])
end

# ╔═╡ 74695bd0-e2b4-4ffd-b089-fdbf81c0b762
md"""
#### Solving Systems of Equations with Augmented Matrices

For this assignment, you will use an augmented matrix to solve a system of linear equations. Refer to [the pdf copy of the Pluto notebook file](https://smccd.instructure.com/courses/54002/files?preview=10121952) we went over in [class on Thursday, Sept. 28](https://smccd.instructure.com/courses/54002/pages/notes-for-09-slash-28-slash-23?module_item_id=3639814), where I showed you how it's done.

Be sure that your augmented matrix has rational number entries. Use row reduction matrices to get the augmented matrix into reduced row echelon form (the coefficient matrix part of the augmented matrix should reduce to the identity matrix).

After you get a solution to the system, use Julia to verify it as shown in that same Pluto notebook file referenced earlier.
"""

# ╔═╡ c297c1fd-6ada-4eb5-a618-d565f7ac5a1a
md"""
###### The function *equations(n)*

I created a function *equations(n)* which returns a system of n equations in n unknowns ``x_1, x_2,\ldots, x_n``.

For this assignment, you'll be solving a system of three linear equations, so enter
```julia
equations(3)
```
in a cell and evaluate it to get your system of three equations. After you've created your augmented matrix of rational values, be sure to disable the cell containing the function. That way it becomes frozen and can't be accidentally changed.

**Extra Credit:** You can earn extra credit points by solving a system of four equations.
"""

# ╔═╡ 8da9ac7e-5fc7-4eaa-82cb-3f403a55da88
equations(4)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
LaTeXStrings = "~1.3.0"
Latexify = "~0.16.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "d5ed385503660414a3803ae7f33df86c624af27f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─2fa9b038-5e2c-11ee-2a7f-29e9335e0135
# ╟─076ea2b4-fcc0-4dbe-8242-ef7344e982f0
# ╟─3e5e4708-7e41-4d62-8f50-806f8d47fe67
# ╟─74695bd0-e2b4-4ffd-b089-fdbf81c0b762
# ╟─c297c1fd-6ada-4eb5-a618-d565f7ac5a1a
# ╠═8da9ac7e-5fc7-4eaa-82cb-3f403a55da88
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
