

"""
This function makes an elliptic solenoid and calculates the inductance.
    The function calculates the flux through the center of each loop's cross-section
        Lₛ = Φ/I (self inductance = flux per amp)
        Lₘ₂₁ = (Φ₂₁ ⋅ Φ₂₂/|Φ₂₂|) /I
            -mutual inductance from the 1ˢᵗ wire on 2ⁿᵈ is flux from 1 in the 2ⁿᵈ's cross section
            that is aligned with the flux from the second on itself. All per ampere
        Lₘ₂₁ = Lₘ₁₂
        Total L = ∑ᵢ∑ⱼ(Lᵢⱼ)
    The inputs are:
    N - Number of turns
    r₁ - one radius of the ellipse
    r₂ - Second radius of ellipse
    L - length of solenoid (axially)
    kwargs:
    MeasLayers - Number of concentric layers to calc field
                    The field is calculated over concentric ellipses (linearly spaced)
                    This is the number of concentric test point layers
    NLayers -  Number of test points per layer
    NPtsPath  - number of points that each wire loop will be discretized into,
    PlotOn deternimes if the coil will be plotted as it builds
        
"""
function eval_Solenoid_Induct(
    N, #Number of turns
    r₁, #Radius 1 of ellipse
    r₂,#Radius 2 of ellipse
    L; #Length of solenoid in meters
    DownSampleFac=1,PlotOn=false,NPtsPath=100,NLayers=20)

    if PlotOn
        pygui(true)
        gcf()
    end
    PointPath_SingleLoop =  MakeEllip(r₁, r₂;NPts=NPtsPath)
    
    LMat = zeros(N, N)
    if PlotOn
        scatter3D(PointPath_SingleLoop[:,1], PointPath_SingleLoop[:,2], PointPath_SingleLoop[:,3])
    end
    NewPP = MakeEllip(r₁,r₂;Center = [0, 0, 1/ N * L],NPts = NPtsPath)

    if PlotOn
        scatter3D(NewPP[:,1], NewPP[:,2], NewPP[:,3])
    end

    Mut, Self, SavedBSelfArr = Mutual_L_TwoLoops(PointPath_SingleLoop, NewPP;DownSampleFac=DownSampleFac,MeasLayers=NLayers,MinThreshold=1e-10,IncludeWireInduct=false,SaveΦ₁=true)
    LMat .+= Self .* OffDiagOnes(N, 1)
    LMat .+= Mut .* OffDiagOnes(N, 2)
    for i in 3:N
        NewPP = MakeEllip(r₁,r₂;Center = [0, 0, (i-1)/ N * L],NPts = NPtsPath)

        if PlotOn
            scatter3D(NewPP[:,1], NewPP[:,2], NewPP[:,3])
        end
        Mut = Mutual_L_TwoLoops(PointPath_SingleLoop, NewPP, SavedBSelfArr;DownSampleFac=DownSampleFac,MeasLayers=NLayers,MinThreshold=1e-10)
        println(Mut)
        
        LMat .+= Mut .* OffDiagOnes(N, i)
        
        println(i)
    end
    LMat[N,N] = LMat[1,1]
    return LMat, sum(LMat)


end

"""
This function makes an elliptic solenoid and calculates the field at the center of the solenoid
        Lₘ₂₁ = Lₘ₁₂
        Total L = ∑ᵢ∑ⱼ(Lᵢⱼ)
    The inputs are:
    N - Number of turns
    r₁ - one radius of the ellipse
    r₂ - Second radius of ellipse
    L - length of solenoid (axially)
    kwargs:

    NPts_Coil  - number of points that each wire loop will be discretized into,
"""
function EvalField_Centered(
    N, #Number of turns
    r₁, #Radius 1 of ellipse
    r₂,#Radius 2 of ellipse
    L; #Length of solenoid in meters
    NPts_Coil=100
    )



            FieldCentered = [
                BiotSav(
                    MakeEllip(
                                r₁, #Radius 1 of ellipse
                                r₂;#Radius 2 of ellipse
                                Center = [0, 0, (Wᵢ - 1) / N * L], #The solenoid's axis is in Z
                                NPts = NPts_Coil, #How discretized the windings are
                                ),
                        [0, 0, L/2];
                        MinThreshold = 0.001,
                        )[3]
                        for Wᵢ in 1:N
                        ]

    return sum(FieldCentered)
end



