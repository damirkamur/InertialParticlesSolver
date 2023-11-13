module IntertialParticlesSolver
    include("utils.jl"); using .Utils
    include("TaskBase.jl"); using .TaskBase
    include("Parsing.jl"); using .Parsing
    include("TaskCylinder.jl"); using .TaskCylinder

    export CylTask, 
        Space,
        Particle, 
        Fiber,
        Grid,
        Point2D,
        Cell2D,
        solve,
        task_from_parse,
        draw,
        Parser,
        get_value,
        draw_field,
        TaskInfo,
        show_info,
        get_xperiod_catch
end
