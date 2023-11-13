module TaskCylinder
    using IntertialParticlesSolver.Utils
    using IntertialParticlesSolver.TaskBase
    using IntertialParticlesSolver.Parsing
    
    import DifferentialEquations as DE
    import Luxor as LX
    import Dates
    
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
    

    mutable struct CylTask{TGrid}
        space::Space
        particles::Vector{Particle}
        N::Int64
        terminated::Int64
        tend::Float64
        max_period::Int64
        St::Float64
        delta::Float64

        grid::TGrid
        periodic::Bool
        fibers::Vector{Fiber}

        _current_fiber::Int64

        name::String
        
        function CylTask(
            space::Space,
            particles::Vector{Particle},
            tend::Real,
            max_period::Int64,
            St::Real,
            delta::Real,
            grid::Union{Grid, Nothing},
            fibers::Vector{Fiber},
            name::String=""
            )
            tend > 0 || error("Время симуляции должно быть больше 0")
            St > 0 || error("Число Стокса должно быть больше 0")
            delta >= 0 || error("Число `delta` должно быть больше или равно 0")
            max_period >= 1 || error("`max_period` должно быть больше или равно 1")
            new{typeof(grid)}(
                space, 
                particles, 
                length(particles), 
                0, 
                float(tend), 
                max_period, 
                float(St), 
                float(delta), 
                grid, 
                grid isa Grid{:Periodic} ? true : false, 
                fibers, 
                0, 
                name
            )
        end
    end

    # x = u[1]; y = u[2]; vx = u[3]; vy = u[4]
    function equations!(du, u, p, t)
        ux, uy = get_value(p[2].grid, u[1], u[2]) 
        du[3] = (ux - u[3]) / p[1]
        du[4] = (uy - u[4]) / p[1] # + vs / p[1]
        du[1] = u[3] 
        du[2] = u[4] 
    end

    function np_condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
        k = 1.0
        for f in integrator.p[2].fibers
            k *= (u[1] - f.x) * (u[1] - f.x) + (u[2] - f.y) * (u[2] - f.y) - (f.r + integrator.p[2].delta) * (f.r + integrator.p[2].delta)
        end

        lx = integrator.p[2].space.x0
        rx = lx + integrator.p[2].space.L
        ly = integrator.p[2].space.y0
        uy = ly + integrator.p[2].space.H
        k *= rx - u[1] 
        k *= u[1] - lx 
        k *= uy - u[2]
        k *= u[2] - ly
        return k
    end

    function p_condition(u, t, integrator)
        k = 1.0
        # Расстояния до волокон
        for f in integrator.p[2].fibers
            x, y, _, _ = _define_conforming_point(integrator.p[5], integrator.p[6], integrator.p[7], integrator.p[8], u[1], u[2])
            fx, fy, _, _ = _define_conforming_point(integrator.p[5], integrator.p[6], integrator.p[7], integrator.p[8], f.x, f.y)
            k *= (x - fx) * (x - fx) + (y - fy) * (y - fy) - (f.r + integrator.p[2].delta) * (f.r + integrator.p[2].delta)
        end 

        # Расстояния до границы расчета
        lx = integrator.p[2].space.x0
        rx = lx + integrator.p[2].space.L * integrator.p[2].max_period
        ly = integrator.p[2].space.y0
        uy = ly + integrator.p[2].space.H * integrator.p[2].max_period
        lx -= integrator.p[2].space.L * (integrator.p[2].max_period - 1)
        ly -= integrator.p[2].space.H * (integrator.p[2].max_period - 1)
        k *= rx - u[1] 
        k *= u[1] - lx 
        k *= uy - u[2]
        k *= u[2] - ly

        return k
    end

    function np_affect!(integrator)
        for f in integrator.p[2].fibers
            if (integrator.u[1] - f.x) * (integrator.u[1] - f.x) + (integrator.u[2] - f.y) * (integrator.u[2] - f.y) - (f.r + integrator.p[2].delta) * (f.r + integrator.p[2].delta) < 1.0e-13
                integrator.p[2].particles[integrator.p[3][1]].terminated = true
                break
            end
        end
        
        DE.terminate!(integrator)
    end

    function p_affect!(integrator)
        for f in integrator.p[2].fibers
            x, y, _, _ = _define_conforming_point(integrator.p[5], integrator.p[6], integrator.p[7], integrator.p[8], integrator.u[1], integrator.u[2])
            fx, fy, _, _ = _define_conforming_point(integrator.p[5], integrator.p[6], integrator.p[7], integrator.p[8], f.x, f.y)
            if (x - fx) * (x - fx) + (y - fy) * (y - fy) - (f.r + integrator.p[2].delta) * (f.r + integrator.p[2].delta) < 1.0e-13
                integrator.p[2].particles[integrator.p[3][1]].terminated = true
                break
            end
        end 
        
        DE.terminate!(integrator)
    end

    
    function solve(t::CylTask)
        if isnothing(t.grid)
            @info "Сетка задачи $(t.name) не задана. Расчет невозможен"
            return 
        end
        tspan = (0.0, t.tend)
        periodic = t.periodic 
        gminx = t.space.x0 
        gminy = t.space.y0 
        gmaxx = gminx + t.space.L 
        gmaxy = gminy + t.space.H

        # iter = 0
        Threads.@threads for i in 1:length(t.particles)
            # iter += 1
            # @info iter
            # t1 = time()
            particle = t.particles[i]
            p = (t.St, t, [0], periodic, gminx, gmaxx, gminy, gmaxy)
            p[3][1] = i
            u0 = [particle.x0, particle.y0, particle.vx0, particle.vy0]
            prob = DE.ODEProblem(equations!, u0, tspan, p)
            condition = periodic ? p_condition : np_condition
            affect! = periodic ? p_affect! : np_affect!
            cb = DE.ContinuousCallback(condition, affect!)
            dt_save = 0.01
            # t2 = time()
            # @info "Подготовка к решению: $(t2 - t1) сек"
            # t1 = time()
            sol = DE.solve(prob, DE.Tsit5(), callback = cb, saveat=tspan[1]:dt_save:tspan[end]+dt_save)
            # t2 = time()
            # @info "Решение: $(t2 - t1) сек"
            # t1 = time()
            particle.ts = sol.t
            particle.xs = [info[1] for info in sol.u]
            particle.ys = [info[2] for info in sol.u]
            # t2 = time()
            # @info "Запись в массивы: $(t2 - t1) сек"
        end
        nothing
    end

    mutable struct TaskInfo
        name::String
        N::Int64
        terminated::Int64
        periodic_terminated::Dict{Tuple, Int64}
        St::Real

        function TaskInfo(task::CylTask)
            info = _collect_info!(task)
            new(info[1], info[2], info[3], info[4], info[5])
        end
    end

    function _collect_info!(t::CylTask)
        name = t.name
        N = t.N
        terminated = sum([p.terminated for p in t.particles])
        periodic_terminated = Dict{Tuple, Int64}()
        for (i, p) in enumerate(t.particles)
            if !p.terminated
                continue
            end
            p.death_period_x, p.death_period_y = _define_period(p.xs[end], p.ys[end], t.space.x0, t.space.L, t.space.y0, t.space.H)
            key = (p.death_period_x, p.death_period_y)
            if haskey(periodic_terminated, key)
                periodic_terminated[key] += 1
            else
                periodic_terminated[key] = 1
            end
        end
        St = t.St
        return (name, N, terminated, periodic_terminated, St)
    end

    function show_info(ti::TaskInfo)
        @info "="^20
        @info "task_name = $(ti.name)"
        @info "N = $(ti.N)"
        @info "terminated = $(ti.terminated)"
        @info "St = $(ti.St)"
        @info "periodic_terminated:"
        for (key, value) in ti.periodic_terminated
            @info "\t$key → $value"
        end
        @info "="^20
        nothing       
    end

    function get_xperiod_catch(ti::TaskInfo, period::Int64)
        catched = 0
        for (key, value) in ti.periodic_terminated
            if key[1] <= period
                catched += value
            end
        end
        return catched
    end
   

    function task_from_parse(
            p::Parser, 
            particles::Vector{Particle}; 
            tend::Real = 10.0, 
            max_period::Int64=5, 
            St::Real = 0.01, 
            N::Int64 = 20, 
            M::Int64 = 20, 
            delta::Real=0, 
            periodic::Bool=false, 
            name::String="")::CylTask
        name = isempty(name) ? p.file_name : name
        L = p.ft.corners[2].x - p.ft.corners[1].x
        H = p.ft.corners[4].y - p.ft.corners[1].y
        space = Space(p.ft.corners[1].x,  p.ft.corners[1].y, L, H)
        grid = if periodic 
            Grid{:Periodic}(p.ftg.points, p.ftg.cells, p.ftg.uxs, p.ftg.uys, N, M)
        else
            Grid{:NonPeriodic}(p.ftg.points, p.ftg.cells, p.ftg.uxs, p.ftg.uys, N, M)
        end
        fibers = p.ft.fibers
        return CylTask(space, particles, tend, max_period, St, delta, grid, fibers, name)
    end

    function draw(
        t, #TODO CylTask 
        desired_size::Int; 
        line_size::Real=1, 
        file_name::String="",
        draw_grid::Bool=true,
        mark_main_cell::Bool=true,
        dx_left::Int64=-1,
        dx_right::Int64=-1,
        dy_left::Int64=-1,
        dy_right::Int64=-1,
        colorize_terminated::Bool=false,
        dark_theme::Bool=false
    )
        maximumx = typemin(Int64)
        maximumy = typemin(Int64)
        minimumx = typemax(Int64)
        minimumy = typemax(Int64)
        for p in t.particles
            maximumx = max(maximumx, max(p.xs...))
            maximumy = max(maximumy, max(p.ys...))
            minimumx = min(minimumx, min(p.xs...))       
            minimumy = min(minimumy, min(p.ys...))
        end
        gminx = t.space.x0 
        gminy = t.space.y0 
        gmaxx = gminx + t.space.L 
        gmaxy = gminy + t.space.H
        realminx = min(minimumx, gminx)
        realmaxx = max(maximumx, gmaxx)
        realminy = min(minimumy, gminy)
        realmaxy = max(maximumy, gmaxy)
        L = t.space.L
        H = t.space.H
        dx_left = dx_left == -1 ? Int(div(gminx - realminx, L)) + 1 : dx_left
        dx_right = dx_right == -1 ? Int(div(realmaxx - gmaxx, L)) + 1 : dx_right
        dy_left = dy_left == -1 ? Int(div(gminy - realminy, H)) + 1 : dy_left
        dy_right = dy_right == -1 ? Int(div(realmaxy - gmaxy, H)) + 1 : dy_right
        Nx = dx_left + dx_right + 1
        Ny = dy_left + dy_right + 1
        # size_coeff
        sc = desired_size / max(Nx * L, Ny * H)
        image_size = (Int(ceil(sc * Nx * L)), Int(ceil(sc * Ny * H))) 
        if isempty(file_name)
            file_name = "src/images/$(string(Dates.format(Dates.now(), "yyyy-mm-dd:H:M:S"))) St=$(t.St) N=$(t.N) terminated=$(t.terminated).png"
        end
        LX.Drawing(image_size[1], image_size[2], file_name)
        LX.setline(line_size)
        LX.background(dark_theme ? "black" : "white")
    
        # Отрисовка сетки
        if draw_grid
            LX.sethue(dark_theme ? "#888888" : "#777777")
            for i = 1:Ny-1
                LX.move(LX.Point(0, i*H*sc))
                LX.line(LX.Point(image_size[1], i*H*sc))
                LX.strokepath()
            end
            for j = 1:Nx-1
                LX.move(LX.Point(L*j*sc, image_size[2]))
                LX.line(LX.Point(L*j*sc, 0))
                LX.strokepath()
            end
        end
        # Основная область
        if mark_main_cell
            LX.sethue(dark_theme ? "#eeeeee" : "#111111")
            LX.setline(line_size+2)
            LX.rect(LX.Point(dx_left*L*sc, dy_right*H*sc), L*sc, H*sc, action = :stroke)
            LX.setline(line_size)
        end
        # Отрисовка волокон
        LX.sethue(dark_theme ? "#ffffff" : "#000000")
        for i = 1:Ny
            for j = 1:Nx
                mirror_x = iseven(j + dx_left)
                mirror_y = iseven(i + dy_left)
                l_dx = (j - 1) * L 
                l_dy = (i - 1) * H
                for f in t.fibers
                    x = (l_dx + (mirror_x ? gmaxx - f.x : f.x - gminx)) * sc
                    y = image_size[2] - (l_dy + (mirror_y ? gmaxy - f.y : f.y - gminy)) * sc
                    r = f.r * sc
                    LX.circle(LX.Point(x, y), r, :fill)
                end 
            end
        end
    
        # Отрисовка частиц
        dx = dx_left*L-gminx
        dy = dy_left*H-gminy
        
        for particle in t.particles  
            LX.sethue(dark_theme ? "#ffffff" : "#000000")
            if colorize_terminated && particle.terminated
                LX.sethue("#ff0000")  
            end 
            LX.move(LX.Point((particle.xs[1] + dx)*sc, image_size[2]-(particle.ys[1] + dy)* sc))
            for i=2:length(particle.xs)
                LX.line(LX.Point((particle.xs[i] + dx)*sc, image_size[2]-(particle.ys[i] + dy)* sc))
            end
            LX.strokepath()
        end
    
        LX.finish()
        nothing
    end

    function draw_field(
        t::CylTask,
        desired_size::Int64=1000;
        dx_left::Int64=0,
        dx_right::Int64=0,
        dy_left::Int64=0,
        dy_right::Int64=0,
        draw_grid::Bool=true,
        mark_main_cell::Bool=true,
        file_name::String="",
        line_size::Real=1.0,
        N::Int64=10,
        point_size::Real=1.0,
        velocity_scale::Real=1.0,
        arrowheadlength::Real=5,
        arrowwidth::Real=2
    )
        if isnothing(t.grid)
            @info "Сетка задачи $(t.name) не задана. Построение поля невозможно"
            return 
        end
        gminx = t.space.x0 
        gminy = t.space.y0 
        gmaxx = gminx + t.space.L 
        gmaxy = gminy + t.space.H
    
        L = t.space.L 
        H = t.space.H
    
        Nx = dx_left + dx_right + 1
        Ny = dy_left + dy_right + 1
        # size_coeff
        sc = desired_size / max(Nx * L, Ny * H)
        image_size = (Int(ceil(sc * Nx * L)), Int(ceil(sc * Ny * H))) 
        if isempty(file_name)
            file_name = "src/images/$(string(Dates.format(Dates.now(), "yyyy-mm-dd:H:M:S"))) St=$(t.St) N=$(t.N) terminated=$(t.terminated).png"
        end
        LX.Drawing(image_size[1], image_size[2], file_name)
        LX.setline(line_size)
        LX.background("white")
    
        # Отрисовка сетки
        if draw_grid
            LX.sethue("#777777")
            for i = 1:Ny-1
                LX.move(LX.Point(0, i*H*sc))
                LX.line(LX.Point(image_size[1], i*H*sc))
                LX.strokepath()
            end
            for j = 1:Nx-1
                LX.move(LX.Point(L*j*sc, image_size[2]))
                LX.line(LX.Point(L*j*sc, 0))
                LX.strokepath()
            end
        end
        # Основная область
        if mark_main_cell
            LX.sethue("#111111")
            LX.setline(line_size+2)
            LX.rect(LX.Point(dx_left*L*sc, dy_right*H*sc), L*sc, H*sc, action = :stroke)
            LX.setline(line_size)
        end
        # Отрисовка волокон
        LX.sethue("#000000")
        for i = 1:Ny
            for j = 1:Nx
                mirror_x = iseven(j + dx_left)
                mirror_y = iseven(i + dy_left)
                l_dx = (j - 1) * L 
                l_dy = (i - 1) * H
                for f in t.fibers
                    x = (l_dx + (mirror_x ? gmaxx - f.x : f.x - gminx)) * sc
                    y = image_size[2] - (l_dy + (mirror_y ? gmaxy - f.y : f.y - gminy)) * sc
                    r = f.r * sc
                    LX.circle(LX.Point(x, y), r, :fill)
                end 
            end
        end
        
        LX.sethue("#ff0000")
        NNx, NNy = L * Nx > H * Ny ? (N, Int(floor(N * H / L))) : (Int(floor(N * L / H)), N)
        NNx = max(NNx, 2)
        NNy = max(NNy, 2)
        dx = dx_left*L
        dy = dy_left*H
        for i=1:Ny
            for j=1:Nx
                for ly in LinRange(0, H, NNy)
                    for lx in LinRange(0, L, NNx)
                        x = lx + gminx + (j-1)*L-dx
                        y = ly + gminy + (i-1)*H-dy
                        ix = (x - gminx + dx) * sc
                        iy = image_size[2] - (y - gminy + dy) * sc
                        vel = get_value(t.grid, x, y)
                        x1 = ((vel[1] * velocity_scale + x) - gminx+dx) * sc
                        y1 = image_size[2] - ((vel[2] * velocity_scale + y) - gminy+dy) * sc
                        if ix != x1 || iy != y1
                            LX.arrow(LX.Point(ix, iy), LX.Point(x1, y1), arrowheadlength=arrowheadlength, linewidth=arrowwidth)
                        end
                        LX.circle(LX.Point(ix, iy), point_size, :fill)
                    end
                end
            end
        end    
        LX.finish()
        nothing
    end
end