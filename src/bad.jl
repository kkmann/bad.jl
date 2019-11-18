module bad

import Base.show

include("util.jl")
export Futility, Efficacy

include("Design.jl")
export Design, n, n1, c

end # module
