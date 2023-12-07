### A Pluto.jl notebook ###
# v0.19.29

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

# ╔═╡ 53cd1e03-7eeb-4046-99b0-ae8aee0c9de6
md"""
##### First I will do a system of 3 equations
"""

# ╔═╡ b557f6ef-48fc-4fd0-811a-248f682becf0
A_aug = [-8//1 -2//1 5//1 6//1; 7//1 2//1 -9//1 5//1; -8//1 9//1 1//1 9//1]

# ╔═╡ 2b186ed4-209b-4ea6-ae7b-41f680568b52
latexify(A_aug)

# ╔═╡ a7900971-e1d3-4959-bffd-c8a15c967f49
E1 = 
[
1//1 0//1 0//1;
0//1 1//1 0//1;
-1//1 0//1 1//1
]

# ╔═╡ 7357ca5d-f471-4c21-bce7-29e92c6421c8
E2 = 
[
1//1 0//1 0//1;
7//8 1//1 0//1;
0//1 0//1 1//1
]

# ╔═╡ dd70242e-243b-422d-8f06-18f16af9f271
E3 = 
[
1//1 0//1 0//1;
0//1 1//1 0//1;
0//1 -44//1 1//1
]

# ╔═╡ fd9f08f3-6860-4bc0-8cb5-8d9cc68392f7
E4 = 
[
1//1 0//1 0//1;
0//1 1//1 74//3192;
0//1 0//1 1//1
]

# ╔═╡ fe0df437-7da0-48a5-917a-ed6d707a0607
E5 = 
[
1//1 0//1 -10//399;
0//1 1//1 0//1;
0//1 0//1 1//1
]

# ╔═╡ e3c0a9b2-8bbb-48ba-a4b3-001419ca647d
E6 = 
[
1//1 8//1 0//1;
0//1 1//1 0//1;
0//1 0//1 1//1
]

# ╔═╡ 4e3505d6-588f-482b-908a-781a502670d1
B =
[
-1//8 0//1 0//1;
0//1 4//1 0//1;
0//1 0//1 2//399
]

# ╔═╡ 33f96b9d-06a5-4eb5-ad7b-81e001abea4d
E = E6*E5*E4*E3*E2*E1

# ╔═╡ 8980c1ac-1104-4c67-8957-0fdfb5467c7a
A_aug_rref = B*E*A_aug

# ╔═╡ fd4029e0-d664-4402-bce5-20f5c9559fbf
latexify(A_aug_rref)

# ╔═╡ b34756d6-0adb-44c6-8ec3-d04f48963072
md"""
$x_1 = \frac{-115}{57}$
$x_2 = \frac{-31}{57}$
$x_3 = \frac{-128}{57}$
"""

# ╔═╡ 01f59290-a096-4fa0-8d51-721e158ecb05
md"""
###### Now to verify the solution:
"""

# ╔═╡ b0da851e-1d4f-4013-ada9-6a4452785fe4
latexify((-115//57)*A_aug[:,1] + (-31//57)*A_aug[:,2] + (-128//57)*A_aug[:,3])

# ╔═╡ 6a058a39-6d81-40ad-9e98-fa375d1a316d
md"""
As you can see, the result matches the $\vec{x}$ coefficient vector from the beginning
"""

# ╔═╡ 3d70ad8d-8581-483d-a2ce-cf3755233b10
md"""
##### EC: Now I will do a system of four equations
"""

# ╔═╡ e3f6a882-b75e-4bb5-b630-71bd4cdaecfb
R_aug = 
[
-7//1 -1//1 3//1 -5//1 -3//1;
3//1 -3//1 -5//1 -1//1 7//1;
-1//1 6//1 8//1 -5//1 2//1;
-6//1 -7//1 -9//1 5//1 9//1
]

# ╔═╡ 3e70ac2d-3fa2-4a50-9f53-d564a5c22933
e1 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 2//1 0//1 1//1	
]

# ╔═╡ bd5a2271-7c58-4c47-92c7-99369b3eb763
e2 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 1//3 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ 7b8ecaea-31d0-4b8a-8160-f0868a47d934
e3 = 
[
1//1 0//1 0//1 0//1;
3//7 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ aca7dc20-aa75-4c4e-b958-fcac5abedad9
e4 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 13//5 1//1	
]

# ╔═╡ b1e23d35-c7b6-4e3b-9417-f6cc600e981f
e5 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 35//24 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ fa36418a-6d6f-4b5d-83ff-24bc1ec85162
e6 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 (12*38)//(11*15) 1//1	
]

# ╔═╡ 72beb8b4-830c-486d-b235-300a6c82e614
e7 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 (-11*119)//(421*12);
0//1 0//1 0//1 1//1	
]

# ╔═╡ cb3ea9a1-4ee4-411a-9058-09b9b9cd25cf
e8 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 (-11*22)//(421*7);
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ 3a8d334d-cf1e-4e3d-a646-421b41c20808
e9 = 
[
1//1 0//1 0//1 (-11*5)//(421);
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ a0e6cf78-900e-4b01-b838-c111e6cae43d
e10 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 (12*26)//(11*7) 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ 22e1101b-58fe-413b-ba61-c6d7b9a39439
e11 = 
[
1//1 0//1 (-12*3)//(11) 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ 813a109d-7d14-4a31-9377-b3586ddd12f2
e12 = 
[
1//1 (-7)//(24) 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1	
]

# ╔═╡ 5adb0e45-ed25-47a5-893b-4397d768b2d0
e = e12*e11*e10*e9*e8*e7*e6*e5*e4*e3*e2*e1

# ╔═╡ e31cdeda-4643-4d50-8135-e7fcd0fe95a0
b = 
[
-1//7 0//1 0//1 0//1;
0//1 -7//24 0//1 0//1;
0//1 0//1 12//11 0//1;
0//1 0//1 0//1 -11//421	
]

# ╔═╡ 0ed54c1a-4689-425a-b4c6-999de4ec533d
R_aug_rref = b * e * R_aug

# ╔═╡ 9dce56b9-d1b6-4d08-96c7-17b8c1b60cfa
latexify(R_aug_rref)

# ╔═╡ 05935090-483b-4f71-b132-2c57e8878175
md"""
$x_1 = \frac{-690}{421}$
$x_2 = \frac{2625}{421}$
$x_3 = \frac{-2426}{421}$
$x_4 = \frac{-762}{421}$
"""

# ╔═╡ 492ae0ab-5b4f-498d-b7d5-5527a5e9be34
md"""
###### Now to verify the solution:
"""

# ╔═╡ ee55be3a-6bd8-4ac6-a4cf-e7f5910e9cf3
latexify((-690//421)*R_aug[:,1] + (2625//421)*R_aug[:,2] + (-2426//421)*R_aug[:,3] + (-762//421)*R_aug[:,4])

# ╔═╡ f2c71b9c-5f0a-4a87-a345-d05be22e7347
md"""
As you can see, the solution is correct as this matrix matches the original coefficient matrix
"""

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

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "d5ed385503660414a3803ae7f33df86c624af27f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

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
version = "5.8.0+0"
"""

# ╔═╡ Cell order:
# ╟─2fa9b038-5e2c-11ee-2a7f-29e9335e0135
# ╟─076ea2b4-fcc0-4dbe-8242-ef7344e982f0
# ╟─3e5e4708-7e41-4d62-8f50-806f8d47fe67
# ╟─74695bd0-e2b4-4ffd-b089-fdbf81c0b762
# ╟─c297c1fd-6ada-4eb5-a618-d565f7ac5a1a
# ╠═8da9ac7e-5fc7-4eaa-82cb-3f403a55da88
# ╟─53cd1e03-7eeb-4046-99b0-ae8aee0c9de6
# ╠═b557f6ef-48fc-4fd0-811a-248f682becf0
# ╠═2b186ed4-209b-4ea6-ae7b-41f680568b52
# ╠═a7900971-e1d3-4959-bffd-c8a15c967f49
# ╠═7357ca5d-f471-4c21-bce7-29e92c6421c8
# ╠═dd70242e-243b-422d-8f06-18f16af9f271
# ╠═fd9f08f3-6860-4bc0-8cb5-8d9cc68392f7
# ╠═fe0df437-7da0-48a5-917a-ed6d707a0607
# ╠═e3c0a9b2-8bbb-48ba-a4b3-001419ca647d
# ╠═4e3505d6-588f-482b-908a-781a502670d1
# ╠═33f96b9d-06a5-4eb5-ad7b-81e001abea4d
# ╠═8980c1ac-1104-4c67-8957-0fdfb5467c7a
# ╠═fd4029e0-d664-4402-bce5-20f5c9559fbf
# ╠═b34756d6-0adb-44c6-8ec3-d04f48963072
# ╠═01f59290-a096-4fa0-8d51-721e158ecb05
# ╠═b0da851e-1d4f-4013-ada9-6a4452785fe4
# ╠═6a058a39-6d81-40ad-9e98-fa375d1a316d
# ╟─3d70ad8d-8581-483d-a2ce-cf3755233b10
# ╠═abeb7093-5afc-4320-aade-06cb12c181f5
# ╠═e3f6a882-b75e-4bb5-b630-71bd4cdaecfb
# ╠═3e70ac2d-3fa2-4a50-9f53-d564a5c22933
# ╠═bd5a2271-7c58-4c47-92c7-99369b3eb763
# ╠═7b8ecaea-31d0-4b8a-8160-f0868a47d934
# ╠═aca7dc20-aa75-4c4e-b958-fcac5abedad9
# ╠═b1e23d35-c7b6-4e3b-9417-f6cc600e981f
# ╠═fa36418a-6d6f-4b5d-83ff-24bc1ec85162
# ╠═72beb8b4-830c-486d-b235-300a6c82e614
# ╠═cb3ea9a1-4ee4-411a-9058-09b9b9cd25cf
# ╠═3a8d334d-cf1e-4e3d-a646-421b41c20808
# ╠═a0e6cf78-900e-4b01-b838-c111e6cae43d
# ╠═22e1101b-58fe-413b-ba61-c6d7b9a39439
# ╠═813a109d-7d14-4a31-9377-b3586ddd12f2
# ╠═5adb0e45-ed25-47a5-893b-4397d768b2d0
# ╠═e31cdeda-4643-4d50-8135-e7fcd0fe95a0
# ╠═0ed54c1a-4689-425a-b4c6-999de4ec533d
# ╠═9dce56b9-d1b6-4d08-96c7-17b8c1b60cfa
# ╠═05935090-483b-4f71-b132-2c57e8878175
# ╠═492ae0ab-5b4f-498d-b7d5-5527a5e9be34
# ╠═ee55be3a-6bd8-4ac6-a4cf-e7f5910e9cf3
# ╠═f2c71b9c-5f0a-4a87-a345-d05be22e7347
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
