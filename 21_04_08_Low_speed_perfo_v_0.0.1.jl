### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 20454a26-4719-4c30-9e46-483c27eb630d
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
			Pkg.PackageSpec(name="PlutoUI", version="0.7"), 
			Pkg.PackageSpec(name="Plots", version="1.10"), 
			Pkg.PackageSpec(name="Colors", version="0.12"),
			Pkg.PackageSpec(name="ColorSchemes", version="3.10"),
			Pkg.PackageSpec(name="LaTeXStrings"),
			])

	using PlutoUI, Plots, Colors, ColorSchemes, LaTeXStrings
	using Statistics, LinearAlgebra  # standard libraries
end

# ╔═╡ e6461173-0482-4912-964d-f349543ed623
PlutoUI.LocalResource("./images/logo.png")

# ╔═╡ f384cb30-50af-4ad5-8986-2e11ed5a5e1d
Markdown.MD(Markdown.Admonition("danger", "DISCLAIMER. STEADY AND LEVEL FLIGHT V0.0.1", [md" This notebook is intended *solely for academic purposes*, It **should not be used** in real operational environments or for aircraft design purposes.  Report issues and find the latest version here  [📡](https://github.com/flt-acdesign/Low_speed_AC_performance)  "]))

# ╔═╡ 6f320fcf-346d-4260-ad63-36269b9de1eb
md"### Set Operating point for calculations ✈  "

# ╔═╡ addf6fa0-3335-4250-9dcb-eb9e3e87b4af
md" Operating Point:       TAS(m/s) =  $(@bind TAS_op NumberField(1:1:340, default=70))  ······   Altitude(m) =  $(@bind Alt_op NumberField(0:100:11000, default=4000))      "

# ╔═╡ d79c73d4-9889-4feb-8eb8-58583dfcc04c
md"### Define aircraft parameters and status   "

# ╔═╡ 22aa1e3d-5265-4e35-90bb-b146954efcf5
md" Wing area in *m^2* (**Sw**) =  $(@bind Sw NumberField(1:1:1000, default=20))     Aircraft maximum lift coefficient (**CLmax**) =  $(@bind CLmax NumberField(0.1:.1:7, default=2))      "

# ╔═╡ de07547d-4c70-4793-96f2-8dfb2379ac54
md" Aircraft mass in *kg* (**M**) =  $(@bind Mass NumberField(0.1:10:600000, default=8000)) ······ **CD0** =  $(@bind CD0 NumberField(0.0:.005:.1, default=.02))   "

# ╔═╡ 98eeefa9-5c98-4c1c-9d00-f494d3eb095a
md"  Aircraft weight in **N** = $(round(Int, Mass*9.81)) ···· Wing Loading (Kgf/m^2) = $(round((Mass/Sw); digits=1))   "

# ╔═╡ 701806df-20fe-4181-b9d0-789f0d4a9944
md" AC/DC (parasit.) = $(round(Int, CD0*10000))  "

# ╔═╡ 508a4cb4-5b8b-40b2-a8bd-33f761a06281
md" e =  $(@bind Oswald NumberField(0.1:.1:1.5, default=.85))  Aspect Ratio =  $(@bind AR NumberField(1:1:30, default=10))  "

# ╔═╡ ae28f744-625b-4fac-a8d2-c74855c752ea
md"### Aircraft Aerodynamic functions"

# ╔═╡ 27007d5a-6f85-4ff4-9185-ab1e0df69eea
md"### ISA+0 Atmosphere functions"

# ╔═╡ 13b762d4-366a-47a5-ad77-a18e2b187b78
begin

# ISA+0 Atmosphere functions
	
	
# NOTE: the text below, with exactly the format used, corresponds to the Julia "docstrings" standard. The documentation needs to be exactly on the line above the function definition. the "Live docs" button at the bottom right of Pluto will show the documentation of the function when the cursor is over the function name anywhere in the code
	
#_________________________________________________________________________________
"""
    ρ(h)

ISA+0 Density in Kg/m^3 as a function of altitude (h) with h in meters (Troposphere). 

Note ρ is written with \\rho<tab>

# Examples
```julia-repl
julia> ρ(0)
1.225
```
"""
ρ(h) = 1.225 * (1-0.0000226 * h) ^4.256

#_________________________________________________________________________________
	
	
#_________________________________________________________________________________

"""
    p(h)

ISA+0 Pressure in Pa as a function of altitude (h) with h in meters (Troposphere).

# Examples
```julia-repl
julia> p(0)
101325
```
"""	
p(h) = 101325 * (1-0.0000226 * h) ^5.256
	
#_________________________________________________________________________________	
	
	
#_________________________________________________________________________________
"""
    T(h)

ISA+0 Temperature in K as a function of altitude (h) with h in meters (Troposphere).

# Examples
```julia-repl
julia> T(0)
288.15
```
"""		
T(h) = 288.15 -6.5 * h /1000
#_________________________________________________________________________________	
	
	
#_________________________________________________________________________________
"""
    a(h)

ISA+0 **speed of sound** in m/s as a function of altitude (**h**) with h in meters (Troposphere). 
	
Note `√` is written with \\sqrt<tab>.

# Examples
```julia-repl
julia> a(0)
340.2626485525556
```
"""			
a(h) = √(1.4*287*T(h))	
#_________________________________________________________________________________
	
		
#_________________________________________________________________________________
"""
    M(TAS, h)

Mach number for a given *True Air Speed (TAS)* in m/s and altitude *h* in meters (troposphere) assuming ISA+0 conditions.

# Examples
```julia-repl
julia> M(70,0)
0.20572343246540056
```
"""			
M(TAS, h) = TAS / a(h)
#_________________________________________________________________________________	
	
	
#_________________________________________________________________________________
"""
    q(TAS, h)

ISA+0 **Dynamic pressure (q)** in Pa as a function of True Air Speed in m/s and altitude (**h**) with h in meters (Troposphere).

# Examples
```julia-repl
julia> q(70,0)
3001.25
```
"""		
q(TAS,h) = .5 * ρ(h)* TAS^2
#_________________________________________________________________________________
	

#_________________________________________________________________________________
"""
    TAS2EAS(v, h)

Convert a True Air Speed value (in any units) to Equivalent Air Speed (in the same units) assuming ISA+0 conditions.

# Examples
```julia-repl
julia> TAS2EAS(100, 3000)
86.12224886708844
```
"""		
TAS2EAS(v, h) = v * (ρ(h)/ρ(0))^.5
#_________________________________________________________________________________	

	
#_________________________________________________________________________________
"""
    EAS2TAS(v, h)

Convert an Equivalent Air Speed value (in any units) to True Air Speed (in the same units) assuming ISA+0 conditions.

# Examples
```julia-repl
julia> EAS2TAS(100, 3000)
116.11401387616912
```
"""		
EAS2TAS(v, h) = v * (ρ(h)/ρ(0))^-.5
#_________________________________________________________________________________	

	
# **** TODO list ****
	
# viscosity
	
# Re/m
	
# Extend to Stratosphere
	
# Add temperature shift
	
	
	
	
	
end	;

# ╔═╡ 19816267-988f-45d5-8c39-bedc11d76e12
md"TAS(m/s) = $(TAS_op) ······ TAS(kt) = $(round(TAS_op*1.94384; digits = 1)) ······ 
EAS(m/s) = $(round(TAS_op*(ρ(Alt_op)/ρ(0))^.5; digits = 1)) ······ EAS(kt) =  $(round(TAS_op*(ρ(Alt_op)/ρ(0))^.5*1.94384; digits = 1))"

# ╔═╡ 29124d03-7ae5-49f1-8aff-272bc9f3d5cd
md"Mach no =  $(round(M(TAS_op, Alt_op); digits = 2))  ······ Altitude (m) =  $(Alt_op) ······  Altitude (ft) =  $(round(Int,Alt_op*3.28084))     "

# ╔═╡ bef4363e-34ef-499b-85cc-eead54a8ede2
begin

# Aircraft aerodynamic coefficients, drag, power required and helper functions
	
	
# NOTE: the text below, with exactly the format used, corresponds to the Julia "docstrings" standard. The documentation needs to be exactly on the line above the function definition. the "Live docs" button at the bottom right of Pluto will show the documentation of the function when the cursor is over the function name anywhere in the code

#_________________________________________________________________________________
"""
    Vs1gTAS(W, h, CLmax, Sw)

Calculate stall speed as TAS at 1g from weight (W) in Newtons, altitude (h) in meters, aircraft maximum lift coefficient (CL) and wing reference area (Sw) in m^2

# Examples
```julia-repl
julia> Vs1gTAS(80000, 4000, 2.1, 20)
48.241248922681834
```
"""
Vs1gTAS(W, h, CLmax, Sw) = ((W)/(ρ(h)*CLmax*Sw))^0.5

#_________________________________________________________________________________

	
#_________________________________________________________________________________
"""
    Vs1gEAS(W, CLmax, Sw)

Calculate stall speed as EAS at 1g from weight (W) in Newtons, aircraft maximum lift coefficient (CL) and wing reference area (Sw) in m^2

# Examples
```julia-repl
julia> Vs1gEAS(80000, 2.1, 20)
39.432317676705956
```
"""
Vs1gEAS(W, CLmax, Sw) = ((W)/(ρ(0)*CLmax*Sw))^0.5

#_________________________________________________________________________________	
	

#_________________________________________________________________________________
"""
    CL(TAS, h, Sw, W) 

Calculate aircraft lift coefficient (CL) from True Air Speed (TAS) in m/s, altitude (h) in m, wing reference area (Sw) in m^2 and aircraft weight (W) in Newtons
	
# Examples
```julia-repl
julia> CL(70, 4000, 70, 80000)
0.569930962682486
```
"""
CL(TAS, h, Sw, W) = W / (.5 * ρ(h) * TAS^2 * Sw)

#_________________________________________________________________________________	
	
	
#_________________________________________________________________________________
"""
    CDi(CL, AR, e)

Calculate aircraft induced drag coefficient (CDi) from aircraft lift coefficient (CL), wing aspect ratio (AR) and Oswald factor (e)
	
# Examples
```julia-repl
julia> CDi(.6, 10, .9)
0.012732395447351627
```
"""
CDi(CL, AR, e) = CL^2 /(π * AR * e)

#_________________________________________________________________________________	
			

#_________________________________________________________________________________
"""
    Drag_induced(TAS, h, Sw, e, W, AR)

Calculate aircraft induced drag in Newtons from true air speed (TAS) in m/s, altitude (h) in meters, wing reference area (Sw) in m^2, Oswald factor (e), aircraft weight (W) in Newtons and wing aspect ratio (AR).
	
# Examples
```julia-repl
julia> Drag_induced(70, 4000, 40, .85, 60000, 10) 
1680.7534663878039
```
"""
Drag_induced(TAS, h, Sw, e, W, AR) = .5 * ρ(h) * TAS^2 * CDi(CL(TAS, h, Sw, W), AR, e) * Sw

#_________________________________________________________________________________		
	
#_________________________________________________________________________________
"""
    Drag_parasitic(TAS, h, CD0, Sw)

Calculate aircraft parasitic drag in Newtons from true air speed (TAS) in m/s, altitude (h) in meters, aircraft zero lift drag coefficient (CD0) and wing reference area (Sw) in m^2.
	
# Examples
```julia-repl
julia> Drag_parasitic(100, 6000, .03, 60)
5929.753076455423
```
"""
Drag_parasitic(TAS, h, CD0, Sw) = .5 * ρ(h) * TAS^2 * CD0 * Sw

#_________________________________________________________________________________		

#_________________________________________________________________________________
"""
    Thrust_required(TAS, h, CD0, Sw, e, W, AR)

Calculate aircraft thrust required for steady and level flight (numerically equal to the total drag = Induced_drag + Parasitic_drag) in Newtons from true air speed (TAS) in m/s, altitude (h) in meters, aircraft zero lift drag coefficient (CD0), wing reference area (Sw) in m^2, Oswald factor (e), aircraft weight (W) in Newtons and wing aspect ratio (AR).
	
# Examples
```julia-repl
julia> Thrust_required(130, 2000, .02, 30, .85, 60000, 10)
5629.534422614251
```
"""
Thrust_required(TAS, h, CD0, Sw, e, W, AR) = Drag_parasitic(TAS, h, CD0, Sw) + Drag_induced(TAS, h, Sw, e, W, AR)

#_________________________________________________________________________________	
	
	
#_________________________________________________________________________________
"""
    Power_required(TAS, h, CD0, Sw, e, W, AR)

Calculate aircraft power required for steady and level flight in W from true air speed (TAS) in m/s, altitude (h) in meters, aircraft zero lift drag coefficient (CD0), wing reference area (Sw) in m^2, Oswald factor (e), aircraft weight (W) in Newtons and wing aspect ratio (AR).
	
# Examples
```julia-repl
julia> Thrust_required(130, 2000, .02, 30, .85, 60000, 10)
731839.4749398526
```
"""
Power_required(TAS, h, CD0, Sw, e, W, AR) = TAS * Thrust_required(TAS, h, CD0, Sw, e, W, AR)
#_________________________________________________________________________________	

	
	
end;

# ╔═╡ da899a54-70a9-4b46-b470-0c9890e5f2de
 Vs1g =   Vs1gTAS(Mass*9.81, Alt_op, CLmax, Sw)

# ╔═╡ c755a5f9-8e92-4612-bfa1-ec543cd66d97
md"EAS(m/s) = $(round(TAS_op*(ρ(Alt_op)/ρ(0))^.5; digits = 1)) ······ EAS(kt) =  $(round(TAS_op*(ρ(Alt_op)/ρ(0))^.5*1.94384; digits = 1)) ······  TASstall(m/s) =  $(round(Vs1g; digits = 1)              )         "

# ╔═╡ ce4bf0a4-97c8-4bf0-9140-d1ff3f05410c
begin

TAS_range = 1:5:360  # Define the range of TAS for the x axis (from 1 to 360 in steps of 5)
h_range = [1:250:10000...] # Define the range of altitudes for the y axis (frfom 1 to 10000 metres in steps of 250m(

	
# Initialize plot	
plot( xticks = 0:50:400, yticks = 0:500:11000, leg=true,
	  grid = (:xy, :olivedrab, :dot, 1, .8) , c= :roma) 

# Draw a contout plot with the dynamic pressure as a function of TAS and Altitude	
plot!(contour(TAS_range, h_range, q, fill = true, c= :coolwarm) )

# Draw a boundary showing the stall speed (TAS) as a function of altitude
plot!(((9.81*Mass)./(ρ.(h_range)*CLmax*Sw)).^.5, h_range, label = "Stall speed", lw= 3)

# Draw a reference line for Mach = 0.5 (below it the flow can be assumed incompressible - although there is no incompressible flow in reality). Below, draw additional Mach boundaries as lines for reference
plot!((0.5.*a.(h_range))  , h_range, label = "M = 0.5", lw= 1)
plot!((0.75.*a.(h_range))  , h_range, label = "M = 0.75", lw= 2)
plot!((1.0.*a.(h_range))  , h_range, label = "M = 1", lw= 3)
	
# Draw a circle showing the operating point under study
scatter!([TAS_op],[Alt_op], label = "Operating Point", ms = 4)	
# Draw a label on the operating point
annotate!([TAS_op]  ,[Alt_op+300], Plots.text("✈", 14, (TAS_op > Vs1g ? :yellow : :red), :left))
# Draw values of the operating point	
annotate!([TAS_op]  ,[Alt_op-300], Plots.text("TAS(kt)= "*string(TAS_op), 8, :yellow, :left))	
annotate!([TAS_op]  ,[Alt_op-650], Plots.text("h(m)= "*string(Alt_op), 8, :yellow, :left))		
annotate!([TAS_op]  ,[Alt_op-1000], Plots.text("q(Pa)= "*string(round(Int,q(TAS_op, Alt_op))), 8, :yellow, :left))			
	

# Draw a circle showing the stall speed at the operating altitude
scatter!([Vs1g],[Alt_op], label = "Stall speed kt(TAS)", ms = 4)	
# Draw a label with the stall speed at this altitude, rounded to 1 decimal place and converted to knots (1m/2 = 1.94384 kt)
annotate!([Vs1g+15]  ,[Alt_op+200], Plots.text(string(round((Vs1g*1.94384); digits=1))*" kt", 8, :orange, :center))
	
# Define plot name and axis labels	
xlabel!("TAS (m/s)")  # Set label for x axis
ylabel!("Altitude (m)")  # Set label for y axis
title!("Dynamic pressure contour with stall and Mach boundaries")
	
plot!(size=(680, 400))	# Update plot attributes

end

# ╔═╡ 887e79e9-c63b-4c52-ade5-44a0bcfdfcf8
begin
	
plot( xticks = 10:20:350, yticks = 0:5000:50000, leg=true, size=(680, 400),grid = (:xy, :olivedrab, :dot, .5, .8)     ) # Initialize plot with some basic parameters
	
	# Plotting the data
v1 = (Vs1g:1:350)   # Define series for x axis (from 0 to 10000 in steps of 500)

	
plot!(v1, Drag_parasitic.(v1, Alt_op, CD0, Sw), label = "D_parasitic (N)", linewidth =3)

plot!(v1, Drag_induced.(v1, Alt_op, Sw, Oswald, Mass*9.81, AR) , label = "D_Induced(N)", linewidth =3)
	
plot!(v1, (Drag_parasitic.(v1, Alt_op, CD0, Sw) + Drag_induced.(v1, Alt_op, Sw, Oswald, Mass*9.81, AR))      , label = "Thrust Required (N)", linewidth =3)
	
		
# Final plot attributes
xlabel!("TAS (m/s)")  # Set label for x axis
ylabel!("Thrust required (N)")  # Set label for y axis (wrt: "with respect to")
title!("Thrust required for steady and level flight at $(Alt_op) m")
	
plot!()  # Update plot with all of the above
	
end

# ╔═╡ a57e579f-5c6e-48f6-a390-2d3b7b816372
begin
	
	plot( xticks = 0:1000:11000, yticks = 0:10:100, leg=true, size=(680, 400),grid = (:xy, :olivedrab, :dot, .5, .8)     ) # Initialize plot with some basic parameters
	
	# Plotting the data
	h1 = (1:500:11000)   # Define series for x axis (from 0 to 10000 in steps of 500)
		
	plot!(h1, p.(h1)./p(0)*100, label = "p(h)/p(0) %", linewidth =3)
	plot!(h1, ρ.(h1)./ρ(0)*100, label = "ρ(h)/ρ(0) %", linewidth =3)
	plot!(h1, T.(h1)./T(0)*100, label = "T(h)/T(0) %", linewidth =2, line = :dashdot)
	
	# Final plot attributes
	xlabel!("h (m)")  # Set label for x axis
	ylabel!("% wrt MSL values")  # Set label for y axis (wrt: "with respect to")
	title!("% (value at altitude / value at MSL) of ISA quantities")
	
	plot!()  # Update plot with all of the above
	
end

# ╔═╡ 5afceffa-6e23-422e-81ee-4aee76899d93
md" 🌍 Graphics showing relative variation with respect to Mean Sea Level (MSL) values of pressure, density and temperature with altitude in ISA+0 conditions in troposphere"

# ╔═╡ 7693366f-2c0c-4be3-be85-e2d7d7591977
md"""
> [***For help with plots follow this link***](http://docs.juliaplots.org/latest/tutorial/) (and then come back to the Pluto notebook)
"""

# ╔═╡ 8b1ccfb8-6b6c-4caa-823f-b16b59eee635
md" The code below this point is to set-up the notebook"

# ╔═╡ caaadf91-6ddd-4933-bf1a-98fb11ab0fec
TableOfContents(aside=true)

# ╔═╡ Cell order:
# ╟─e6461173-0482-4912-964d-f349543ed623
# ╟─f384cb30-50af-4ad5-8986-2e11ed5a5e1d
# ╟─6f320fcf-346d-4260-ad63-36269b9de1eb
# ╟─addf6fa0-3335-4250-9dcb-eb9e3e87b4af
# ╟─19816267-988f-45d5-8c39-bedc11d76e12
# ╟─29124d03-7ae5-49f1-8aff-272bc9f3d5cd
# ╟─c755a5f9-8e92-4612-bfa1-ec543cd66d97
# ╟─ce4bf0a4-97c8-4bf0-9140-d1ff3f05410c
# ╟─d79c73d4-9889-4feb-8eb8-58583dfcc04c
# ╟─22aa1e3d-5265-4e35-90bb-b146954efcf5
# ╟─de07547d-4c70-4793-96f2-8dfb2379ac54
# ╟─98eeefa9-5c98-4c1c-9d00-f494d3eb095a
# ╟─701806df-20fe-4181-b9d0-789f0d4a9944
# ╟─508a4cb4-5b8b-40b2-a8bd-33f761a06281
# ╠═887e79e9-c63b-4c52-ade5-44a0bcfdfcf8
# ╟─ae28f744-625b-4fac-a8d2-c74855c752ea
# ╠═da899a54-70a9-4b46-b470-0c9890e5f2de
# ╠═bef4363e-34ef-499b-85cc-eead54a8ede2
# ╟─27007d5a-6f85-4ff4-9185-ab1e0df69eea
# ╠═13b762d4-366a-47a5-ad77-a18e2b187b78
# ╟─a57e579f-5c6e-48f6-a390-2d3b7b816372
# ╟─5afceffa-6e23-422e-81ee-4aee76899d93
# ╟─7693366f-2c0c-4be3-be85-e2d7d7591977
# ╟─8b1ccfb8-6b6c-4caa-823f-b16b59eee635
# ╟─caaadf91-6ddd-4933-bf1a-98fb11ab0fec
# ╟─20454a26-4719-4c30-9e46-483c27eb630d
