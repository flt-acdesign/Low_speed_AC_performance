### A Pluto.jl notebook ###
# v0.14.0

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

# ╔═╡ caaadf91-6ddd-4933-bf1a-98fb11ab0fec
TableOfContents(aside=true)

# ╔═╡ 6f320fcf-346d-4260-ad63-36269b9de1eb
md"### Set Operating point for calculations   "

# ╔═╡ afa9f7ee-d8e5-426c-9df0-3a05641d8fbb
md" Operating Point:       TAS(m/s) $(@bind TAS_op Slider(1:1:340))  Alt(m) $(@bind Alt_op Slider(0:100:11000)) "

# ╔═╡ d79c73d4-9889-4feb-8eb8-58583dfcc04c
md"### Define aircraft parameters and status   "

# ╔═╡ 22aa1e3d-5265-4e35-90bb-b146954efcf5
md" Wing area in *m^2* (**Sw**) =  $(@bind Sw NumberField(1:10:1000, default=20))     Aircraft maximum lift coefficient (CLmax) =  $(@bind CLmax NumberField(0.1:.1:7, default=2))      "

# ╔═╡ de07547d-4c70-4793-96f2-8dfb2379ac54
md" Aircraft mass in *kg* (**M**) =  $(@bind Mass NumberField(0.1:10:600000, default=8000)) ______ CD0 =  $(@bind CD0 NumberField(0.0:.005:.1, default=.04))   "

# ╔═╡ 98eeefa9-5c98-4c1c-9d00-f494d3eb095a
md"  Aircraft weight in **N** = $(round(Int, Mass*9.81)) _____  AC/DC = $(round(Int, CD0*10000))    "

# ╔═╡ d6f27d04-c041-4eee-a219-574dbed118fc
md"  Wing Loading (Kgf/m^2) = $(round((Mass/Sw); digits=1)  ) " 

# ╔═╡ 508a4cb4-5b8b-40b2-a8bd-33f761a06281
md" e =  $(@bind Oswald NumberField(0.1:.1:1.5, default=.85))  Aspect Ratio =  $(@bind AR NumberField(1:1:30, default=10))  "

# ╔═╡ ae28f744-625b-4fac-a8d2-c74855c752ea
md"### Aircraft Aerodynamic functions"

# ╔═╡ 5ccc2dc6-7a26-4273-9a6d-1cafd1c7accb
CDi(_CL, _AR, _e) = _CL^2 /(π * _AR * _e);

# ╔═╡ 27007d5a-6f85-4ff4-9185-ab1e0df69eea
md"### ISA+0 Atmosphere functions"

# ╔═╡ 2a6af729-1805-4222-b134-592d25ff1aa5
md" ρ(h) :  ISA+0 Density in Kg/m^3 as a function of altitude (h) with h in meters (Troposphere). Note ρ is written with `\rho<tab>`"

# ╔═╡ 13b762d4-366a-47a5-ad77-a18e2b187b78
ρ(h) = 1.225 * (1-0.0000226 * h) ^4.256;

# ╔═╡ 4077c54d-2a6e-4e80-975c-caf3825f1bc1
Drag_parasitic(_TAS, _alt, _CD0, _Sw) = .5 * ρ(_alt) * _TAS^2 * CD0 * _Sw;

# ╔═╡ b0c37277-a103-4eef-8ec7-544a6a68cb94
CL(_TAS, _alt, _Sw, _W) = _W / (.5 * ρ(_alt) * _TAS^2 * _Sw);

# ╔═╡ 4116c5b3-64f8-4f31-8040-cc997ba57243
Drag_induced(_TAS, _alt, _CD0, _Sw, _e, _W, _AR) = 	.5 * ρ(_alt) * _TAS^2 * CDi(CL(_TAS, _alt, _Sw, _W), _AR, _e) * _Sw;

# ╔═╡ 8cb8c86f-854b-4b7d-af43-197527f1df3e
# Calculation of stall speed at the operating altitude	
TAS_Vs1g_at_alt(_W, _h, _CLmax, _Sw) = ((_W)/(ρ(_h)*_CLmax*_Sw))^0.5 ;

# ╔═╡ da899a54-70a9-4b46-b470-0c9890e5f2de
 Vs1g =   TAS_Vs1g_at_alt(Mass*9.81, Alt_op, CLmax, Sw)

# ╔═╡ c755a5f9-8e92-4612-bfa1-ec543cd66d97
md"EAS(m/s) = $(round(TAS_op*(ρ(Alt_op)/ρ(0))^.5; digits = 1)) ___ EAS(kt) =  $(round(TAS_op*(ρ(Alt_op)/ρ(0))^.5*1.94384; digits = 1)) ____  TASstall(m/s) =  $(round(Vs1g; digits = 1))         "

# ╔═╡ 887e79e9-c63b-4c52-ade5-44a0bcfdfcf8
begin
	
plot( xticks = 10:10:350, yticks = 0:5000:50000, leg=true, size=(680, 400),grid = (:xy, :olivedrab, :dot, .5, .8)     ) # Initialize plot with some basic parameters
	
	# Plotting the data
v1 = (Vs1g:1:350)   # Define series for x axis (from 0 to 10000 in steps of 500)

	
plot!(v1, Drag_parasitic.(v1, Alt_op, CD0, Sw), label = "D_parasitic (N)", linewidth =3)

plot!(v1, Drag_induced.(v1, Alt_op, CD0, Sw, Oswald, Mass*9.81, AR) , label = "D_Induced(N)", linewidth =3)
	
plot!(v1, (Drag_parasitic.(v1, Alt_op, CD0, Sw) + Drag_induced.(v1, Alt_op, CD0, Sw, Oswald, Mass*9.81, AR))      , label = "Thrust Required (N)", linewidth =3)
	
		
# Final plot attributes
xlabel!("TAS (m)")  # Set label for x axis
ylabel!("Thrust required")  # Set label for y axis (wrt: "with respect to")
title!("Thrust required for steady and level flight")
	
plot!()  # Update plot with all of the above
	
end

# ╔═╡ 8ce71cb5-4f4e-4e02-9999-48cadc5a5e59
md" p(h) :  ISA+0 Pressure in Pa as a function of altitude (h) with h in meters (Troposphere)"

# ╔═╡ 63865400-7a48-4c70-b7ca-1e2c6f4456cb
p(h) = 101325 * (1-0.0000226 * h) ^5.256;

# ╔═╡ c1af0de6-87f8-4def-b948-e48f1c550310
md" T(h) :  ISA+0 Temperature in K as a function of altitude (h) with h in meters (Troposphere)"

# ╔═╡ 378f5ff0-9854-11eb-03ba-b79cd02a308b
T(h) = 288.15 -6.5 * h /1000;

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
md" > Graphics showing relative variation with respect to Mean Sea Level (MSL) values of pressure, density and temperature with altitude in ISA+0 conditions up to 10000 metres of altitude"

# ╔═╡ 7693366f-2c0c-4be3-be85-e2d7d7591977
md"""
> [***For help on plots follow this link***](http://docs.juliaplots.org/latest/tutorial/) (and then come back to the Pluto notebook)
"""

# ╔═╡ 3f5be5de-30eb-4fb3-b598-41fba6be075a
md" a(h) :  ISA+0 **speed of sound** in m/s as a function of altitude (**h**) with h in meters (Troposphere). Note `√` is written with `\sqrt<tab>`"

# ╔═╡ 418c035e-6342-464c-9f7b-2c47767a1ede
a(h) = √(1.4*287*T(h)); 

# ╔═╡ a50e01cb-8585-4f2d-94d6-1fef73a35660
md" q(TAS, h) :  ISA+0 **Dynamic pressure (q)** in Pa as a function of True Air Speed in m/s and altitude (**h**) with h in meters (Troposphere). "

# ╔═╡ 4d35293d-3313-4870-98bf-6eeb31a388e0
q(TAS,h) = .5 * ρ(h)* TAS^2;

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
annotate!([TAS_op]  ,[Alt_op+200], Plots.text("OP", 9, :yellow, :left))
# Draw values of the operating point	
annotate!([TAS_op]  ,[Alt_op-300], Plots.text("TAS(kt)= "*string(TAS_op), 8, :yellow, :left))	
annotate!([TAS_op]  ,[Alt_op-600], Plots.text("h(m)= "*string(Alt_op), 8, :yellow, :left))		

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

# ╔═╡ b99e3354-b1c9-4ba9-b3d0-d2ea5450c7c5
md" M(TAS,h) :  Mach number for a given *True Air Speed (TAS)* in m/s and altitude *h* in meters (troposphere) assuming ISA+0 conditions"

# ╔═╡ 867bbb25-27f9-46af-8f61-46877f77d8ef
M(TAS, h) = TAS / a(h) ;

# ╔═╡ 5502a5c3-3ea4-452c-9704-e49e6434aa40
md"TAS(m/s) = $(TAS_op) ___ TAS(kt) = $(round(TAS_op*1.94384; digits = 1)) ___ Mach no =  $(round(M(TAS_op, Alt_op); digits = 2))  ___ Altitude (m) =  $(Alt_op)    "

# ╔═╡ 8b1ccfb8-6b6c-4caa-823f-b16b59eee635
md" The code below this point is to set-up the notebook"

# ╔═╡ Cell order:
# ╟─caaadf91-6ddd-4933-bf1a-98fb11ab0fec
# ╟─6f320fcf-346d-4260-ad63-36269b9de1eb
# ╟─afa9f7ee-d8e5-426c-9df0-3a05641d8fbb
# ╟─5502a5c3-3ea4-452c-9704-e49e6434aa40
# ╟─c755a5f9-8e92-4612-bfa1-ec543cd66d97
# ╟─ce4bf0a4-97c8-4bf0-9140-d1ff3f05410c
# ╟─d79c73d4-9889-4feb-8eb8-58583dfcc04c
# ╟─22aa1e3d-5265-4e35-90bb-b146954efcf5
# ╟─de07547d-4c70-4793-96f2-8dfb2379ac54
# ╟─98eeefa9-5c98-4c1c-9d00-f494d3eb095a
# ╟─d6f27d04-c041-4eee-a219-574dbed118fc
# ╟─508a4cb4-5b8b-40b2-a8bd-33f761a06281
# ╟─887e79e9-c63b-4c52-ade5-44a0bcfdfcf8
# ╟─ae28f744-625b-4fac-a8d2-c74855c752ea
# ╠═4077c54d-2a6e-4e80-975c-caf3825f1bc1
# ╠═4116c5b3-64f8-4f31-8040-cc997ba57243
# ╠═b0c37277-a103-4eef-8ec7-544a6a68cb94
# ╠═5ccc2dc6-7a26-4273-9a6d-1cafd1c7accb
# ╠═8cb8c86f-854b-4b7d-af43-197527f1df3e
# ╠═da899a54-70a9-4b46-b470-0c9890e5f2de
# ╟─27007d5a-6f85-4ff4-9185-ab1e0df69eea
# ╟─2a6af729-1805-4222-b134-592d25ff1aa5
# ╠═13b762d4-366a-47a5-ad77-a18e2b187b78
# ╟─8ce71cb5-4f4e-4e02-9999-48cadc5a5e59
# ╠═63865400-7a48-4c70-b7ca-1e2c6f4456cb
# ╟─c1af0de6-87f8-4def-b948-e48f1c550310
# ╠═378f5ff0-9854-11eb-03ba-b79cd02a308b
# ╟─a57e579f-5c6e-48f6-a390-2d3b7b816372
# ╟─5afceffa-6e23-422e-81ee-4aee76899d93
# ╟─7693366f-2c0c-4be3-be85-e2d7d7591977
# ╟─3f5be5de-30eb-4fb3-b598-41fba6be075a
# ╠═418c035e-6342-464c-9f7b-2c47767a1ede
# ╟─a50e01cb-8585-4f2d-94d6-1fef73a35660
# ╠═4d35293d-3313-4870-98bf-6eeb31a388e0
# ╟─b99e3354-b1c9-4ba9-b3d0-d2ea5450c7c5
# ╠═867bbb25-27f9-46af-8f61-46877f77d8ef
# ╟─8b1ccfb8-6b6c-4caa-823f-b16b59eee635
# ╟─20454a26-4719-4c30-9e46-483c27eb630d
