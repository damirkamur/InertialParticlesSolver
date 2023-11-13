module TaskBase
    using IntertialParticlesSolver.Utils
    
    export Point2D,
        cross,
        Cell2D,
        Triagle,
        Fiber,
        FiberTaskParse,
        FiberTaskGridParse,
        Grid,
        get_value,
        _define_conforming_point,
        Particle,
        Space

    mutable struct Particle
        x0::Float64
        y0::Float64
        vx0::Float64
        vy0::Float64
        ts::Vector{Float64}
        xs::Vector{Float64}
        ys::Vector{Float64}
        terminated::Bool
        death_period_x::Int64
        death_period_y::Int64
        
        function Particle(x0::Real, y0::Real, vx0::Real, vy0::Real)
            ts=Vector{Float64}[]
            xs=Vector{Float64}[]
            ys=Vector{Float64}[]
            new(float(x0), float(y0), float(vx0), float(vy0), ts, xs, ys, false, 0, 0)
        end
    end

    mutable struct Space
        x0::Float64
        y0::Float64
        L::Float64
        H::Float64
    
        """
            {x0 y0} – нижний левый угол прямоугольника
            {L, H} – Ширина и высота соответственно 
        """
        function Space(x0::Real, y0::Real, L::Real, H::Real)
            L > 0 || error("Размер области должен иметь положительное значение")
            H > 0 || error("Размер области должен иметь положительное значение")
            new(float(x0), float(y0), float(L), float(H))
        end
    end

    mutable struct Point2D
        x::Float64
        y::Float64
        
        Point2D(x::Real, y::Real) = new(float(x), float(y))
    end
    Base.:+(p1::Point2D, p2::Point2D) = Point2D(p1.x + p2.x, p1.y + p2.y)
    Base.:-(p1::Point2D, p2::Point2D) = Point2D(p1.x - p2.x, p1.y - p2.y)
    cross(p1::Point2D, p2::Point2D) = p1.x * p2.y - p2.x * p1.y

    mutable struct Cell2D
        p1::Int64
        p2::Int64
        p3::Int64

        function Cell2D(n::Int64, p1::Int64, p2::Int64, p3::Int64) 
            n == 3 || error("Поддреживаются только треугольные ячейки. получен размер: $n")
            new(p1, p2, p3)
        end
    end

    mutable struct Triagle
        p1::Point2D
        p2::Point2D
        p3::Point2D
        left_up::Point2D
        bottom_right::Point2D

        function Triagle(p1::Point2D, p2::Point2D, p3::Point2D)
            f2 = p2 - p1
            f3 = p3 - p1
            cross(f2, f3) > 0 || error("Треугольник ($p1, $p2, $p3) некорректен")
            xs = [p.x for p in[p1, p2, p3]]
            ys = [p.y for p in[p1, p2, p3]]
            left_up = Point2D(minimum(xs), maximum(ys))
            bottom_right = Point2D( maximum(xs), minimum(ys))
            new(p1, p2, p3, left_up, bottom_right)
        end
    end

    function point_in_triangle(p0::Point2D, p1::Point2D, p2::Point2D, p3::Point2D)::Bool
        d1 = cross(p1 - p0, p2 - p0)
        d2 = cross(p2 - p0, p3 - p0)
        d3 = cross(p3 - p0, p1 - p0)

        has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0)
        has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0)

        return !(has_neg && has_pos)
    end

    point_in_triangle(p::Point2D, t::Triagle)::Bool = point_in_triangle(p, t.p1, t.p2, t.p3)

    function interpolate_into_triangle(p1::Point2D, p2::Point2D, p3::Point2D, v1::Real, v2::Real, v3::Real, f::Point2D)::Float64
        f1 = p1 - f
        f2 = p2 - f
        f3 = p3 - f
        a = abs(cross(p1 - p2, p1 - p3))
        a1 = abs(cross(f2, f3)) / a
        a2 = abs(cross(f3, f1)) / a
        a3 = abs(cross(f1, f2)) / a
        return v1 * a1 + v2 * a2 + v3 * a3
    end

    interpolate_into_triangle(t::Triagle, v1::Real, v2::Real, v3::Real, f::Point2D) = interpolate_into_triangle(t.p1, t.p2, t.p3, v1, v2, v3, f)

    mutable struct Fiber
        x::Float64
        y::Float64
        r::Float64

        Fiber(x::Real, y::Real, r::Real) = new(float(x), float(y), float(r))
    end

    mutable struct FiberTaskParse
        n::Int64
        corners::Vector{Point2D}
        fibers::Vector{Fiber}

        function FiberTaskParse(
                n::Int64,
                corners::Vector{Point2D},
                fibers::Vector{Fiber}
            )
            n == length(fibers) || error("Заданное число волокон не совпадает с их фактическим числом")
            length(corners) == 4 || error("Поддерживается только прямоугольная область")
            new(n, corners, fibers)
        end
    end

    mutable struct FiberTaskGridParse
        npoints::Int64
        points::Vector{Point2D}
        ncells::Int64
        cells::Vector{Cell2D}
        uxs::Vector{Float64}
        uys::Vector{Float64}

        function FiberTaskGridParse(
            points::Vector{Point2D},
            cells::Vector{Cell2D},
            uxs::Vector{Float64},
            uys::Vector{Float64}
            )
            npoints = length(points)
            ncells = length(cells)
            length(uxs) == length(uys) == npoints || 
                error("Кол-во информации по точкам не совпадает с кол-вом точек. ux - $(length(uxs)), uy - $(length(uys)), npoints - $npoints")
            
            s1 = Set()
            for cell in cells
                push!(s1, cell.p1, cell.p2, cell.p3)
            end

            s2 = Set(1:npoints)
            s1 == s2 || error("Не все точки сетки используются или используются несущеуствующие")

            new(
                npoints,
                points,
                ncells,
                cells,
                uxs,
                uys
            )
        end
    end

    mutable struct Grid{Periodic}
        npoints::Int64
        points::Vector{Point2D}

        ncells::Int64
        cells::Vector{Cell2D}
        
        triangles::Vector{Triagle}
        
        hashtable::Dict{Int64, Set{Int64}}
        dx::Float64
        dy::Float64

        uxs::Vector{Float64}
        uys::Vector{Float64}
        
        N::Int64
        M::Int64
        
        minx::Float64
        maxx::Float64
        miny::Float64
        maxy::Float64

        periodic::Bool

        # N, M - кол-во элементов для хэш таблицы
        function Grid{T}(points::Vector{Point2D}, cells::Vector{Cell2D}, uxs::Vector{Float64}, uys::Vector{Float64}, N::Int64, M::Int64) where T
            T in (:Periodic, :NonPeriodic) || error("Типизация Grid :Periodic или :NonPeriodic")
            periodic = (T == :Periodic)
            npoints = length(points)
            ncells = length(cells)
            
            triangles = [Triagle(points[cell.p1], points[cell.p2], points[cell.p3]) for cell in cells]
            minx, maxx, miny, maxy = _define_limits(points)
            hashtable, dx, dy = _generate_hash_table(minx, maxx, miny, maxy, N, M, triangles)

            new{T}(npoints, points, ncells, cells, triangles, hashtable, dx, dy, uxs, uys, N, M, minx, maxx, miny, maxy, periodic)
        end 
        
        function Grid(points::Vector{Point2D}, cells::Vector{Cell2D}, uxs::Vector{Float64}, uys::Vector{Float64}, N::Int64, M::Int64)
            return Grid{:NonPeriodic}(points, cells, uxs, uys, N, M)
        end
    end

    function get_value(g::Grid{:NonPeriodic}, point::Point2D)::Tuple{Float64, Float64}
        hashindex = _define_rect(point.x, point.y, g)
        isempty(g.hashtable[hashindex]) && return (0.0, 0.0)
        
        for ti in g.hashtable[hashindex]
            triangle = g.triangles[ti]
            
            if point_in_triangle(point, triangle)
                cell = g.cells[ti]
                ux1 = g.uxs[cell.p1]
                ux2 = g.uxs[cell.p2]
                ux3 = g.uxs[cell.p3]
                uy1 = g.uys[cell.p1]
                uy2 = g.uys[cell.p2]
                uy3 = g.uys[cell.p3]
                ux = interpolate_into_triangle(triangle, ux1, ux2, ux3, point)
                uy = interpolate_into_triangle(triangle, uy1, uy2, uy3, point)
                return ux, uy
            end
        end
        # @warn "Выход точки за границу области: $point"
        return (0.0, 0.0)
    end

    get_value(g::Grid{:NonPeriodic}, x::Real, y::Real)::Tuple{Float64, Float64} = get_value(g, Point2D(x, y))

    function get_value(g::Grid{:Periodic}, point::Point2D)::Tuple{Float64, Float64}
        newx, newy, flagx, flagy = _define_conforming_point(g.minx, g.maxx, g.miny, g.maxy, point.x, point.y)
        newpoint = Point2D(newx, newy)
        hashindex = _define_rect(newpoint.x, newpoint.y, g)
        isempty(g.hashtable[hashindex]) && return (0.0, 0.0)
        
        for ti in g.hashtable[hashindex]
            triangle = g.triangles[ti]
            
            if point_in_triangle(newpoint, triangle)
                cell = g.cells[ti]
                ux1 = g.uxs[cell.p1]
                ux2 = g.uxs[cell.p2]
                ux3 = g.uxs[cell.p3]
                uy1 = g.uys[cell.p1]
                uy2 = g.uys[cell.p2]
                uy3 = g.uys[cell.p3]
                ux = interpolate_into_triangle(triangle, ux1, ux2, ux3, newpoint)
                uy = interpolate_into_triangle(triangle, uy1, uy2, uy3, newpoint)
                return ux, !((flagx && !flagy)|| (!flagx && flagy)) ? uy : - uy
            end
        end
        # @warn "Выход точки за границу области: $point"
        return (0.0, 0.0)
    end

    get_value(g::Grid{:Periodic}, x::Real, y::Real)::Tuple{Float64, Float64} = get_value(g, Point2D(x, y))

    function _define_limits(points::Vector{Point2D})::Tuple{Float64, Float64, Float64, Float64}
        xs = [p.x for p in points]
        ys = [p.y for p in points]
        minx = minimum(xs)
        maxx = maximum(xs)
        miny = minimum(ys)
        maxy = maximum(ys)
        return minx, maxx, miny, maxy
    end

    function _generate_hash_table(
        minx::Float64, 
        maxx::Float64, 
        miny::Float64, 
        maxy::Float64,
        N::Int64,
        M::Int64,
        triangles::Vector{Triagle})::Tuple{Dict{Int64, Set{Int64}}, Float64, Float64}

        xs = LinRange(minx, maxx, N+1)
        dx = xs[2] - xs[1]
        ys = LinRange(miny, maxy, M+1)
        dy = ys[2] - ys[1]
        
        hashtable = Dict{Int64, Set{Int64}}(
            i => Set{Int64}() for i in 1:(N * M)
        )

        Threads.@threads for j = 1:N
            for i = 1:M
                left_up = Point2D(xs[j], ys[i + 1])
                bottom_right = Point2D(xs[j + 1], ys[i])
                for (index, t) in enumerate(triangles)
                    left_up2 = t.left_up
                    bottom_right2 = t.bottom_right
                    if _intersect_rectangles(left_up, bottom_right, left_up2, bottom_right2)
                        push!(hashtable[j + (i - 1)* N], index)
                    end
                end
            end
        end

        return hashtable, dx, dy
    end

    function _intersect_rectangles(LU1::Point2D, BR1::Point2D, LU2::Point2D, BR2::Point2D)
        # return LU1.y < BR2.y || BR1.y > LU2.y || BR1.x < LU2.x || LU1.x > BR2.x
        LU2.y <= BR1.y && return false
        LU2.x >= BR1.x && return false
        BR2.y >= LU1.y && return false
        BR2.x <= LU1.x && return false
        return true
    end

    function _define_rect_zero(x::Float64, y::Float64, dx::Float64, dy::Float64, N::Int64, M::Int64)::Int64
        ix = Int64(floor(x / dx)+1) 
        iy = Int64(floor(y / dy)+1) 
        ix = min(N, ix)
        iy = min(M, iy)
        return ix + (iy - 1) * N
    end

    function _define_rect(x::Float64, y::Float64, g::Grid)::Int64
        dx = g.dx
        dy = g.dy
        N = g.N
        M = g.M
        x -= g.minx
        y -= g.miny
        return _define_rect_zero(x, y, dx, dy, N, M)
    end
end