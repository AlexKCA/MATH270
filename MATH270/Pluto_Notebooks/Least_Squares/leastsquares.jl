### A Pluto.jl notebook ###
# v0.19.30

using Markdown
using InteractiveUtils

# ╔═╡ 83eebda9-e3a9-46c7-a500-7f75c9940790
using LinearAlgebra

# ╔═╡ 34f64474-7dc7-11ee-226e-b7a32880693e
md"""
##### Miles Kent - MATH 270
# Least Squares Assignment

Want the cubic polynomial $y = ax^3 + bx^2 + cx + d$ of best fit for a set of randomly generated $(x, y)$ points

Minimize square error by finding ideal $a, b, c, d$

Submission requires
1. A plot of the data points together with the plot of the best fitting cubic function.
2. The equation of the cubic function rendered in $\LaTeX$.
3. The least square error rendered in $\LaTeX$.
4. Verication of the Pythagorean Theorem for the vectors $p, e, y$ by showing that $p^T p + e^T e \approx y^T y$. You can do this by showing that $|p^T p + e^T e - y^T y| < 10^{-10}$

"""

# ╔═╡ 99d4b653-267b-4bcb-babe-2374a00c8373
md"""
### Generate random $(x, y)$ data points
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "ac1187e548c6ab173ac57d4e72da1620216bce54"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"
"""

# ╔═╡ Cell order:
# ╠═83eebda9-e3a9-46c7-a500-7f75c9940790
# ╟─34f64474-7dc7-11ee-226e-b7a32880693e
# ╟─99d4b653-267b-4bcb-babe-2374a00c8373
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
