module Parsing
    using IntertialParticlesSolver.Utils
    using IntertialParticlesSolver.TaskBase
    
    export Parser,
        Point2D,
        cross,
        Cell2D,
        Triagle,
        Fiber,
        FiberTaskParse,
        FiberTaskGridParse,
        Grid,
        get_value

    mutable struct Parser
        file_name::String
        dat_file::String
        vtk_file::String
        ft::FiberTaskParse
        ftg::FiberTaskGridParse

        function Parser(file_name::String)
            dat_file = "src/tasks/$file_name.dat"
            vtk_file = "src/tasks/$file_name.vtk"
            isfile(dat_file) && isfile(vtk_file) || error("Файлы `$dat_file` и `$vtk_file` не найдены")
            new(
                file_name,
                dat_file, 
                vtk_file, 
                _parse_dat(dat_file),
                _parse_vtk(vtk_file)
                )
        end
    end

    function _parse_dat(file_name::String)::FiberTaskParse
        state = 0

        corners = Point2D[]
        fibers = Fiber[]
        n = 0
        
        open(file_name, "r") do file           
            while !eof(file)
                line = strip(readline(file))
                new_state = _define_dat_state(line)
                state = iszero(new_state) ? state : new_state
                !iszero(new_state) && continue

                state == 3 && continue
                if state == 1
                    try
                        data = split(line, r" +") .|> Meta.parse .|> eval
                        length(data) == 2 || data isa Vector{<:Real} || error("")
                        push!(corners, Point2D(data[1], data[2]))
                    catch
                        error("Ошибка парсинга раздела `! corners (x y)`")
                    end
                elseif state == 2
                    try
                        data = split(line, r" +") .|> Meta.parse .|> eval
                        length(data) == 2 || data isa Vector{Int64} || error("")
                        n = data[1]
                    catch
                        error("Ошибка парсинга раздела `! fibers (n nx)`")
                    end            
                elseif state == 4
                    try
                        data = split(line, r" +") .|> Meta.parse .|> eval
                        length(data) == 3 || data isa Vector{<:Real} || error("")
                        push!(fibers, Fiber(data[1], data[2], data[3]))
                    catch
                        error("Ошибка парсинга раздела `! fibers (x y r)`")
                    end   
                end
            end
        end
     
        return  FiberTaskParse(n, corners, fibers)
    end

    function _define_dat_state(line)::Int64
        contains(line, "! corners (x y)") && return 1
        contains(line, "! fibers (n nx)") && return 2
        contains(line, "! fibers (m1 m2)") && return 3
        contains(line, "! fibers (x y r)") && return 4
        return 0            
    end

    function _parse_vtk(file_name::String)::FiberTaskGridParse
        npd = 0

        points = Point2D[]
        cells = Cell2D[]
        uxs = Float64[]
        uys = Float64[]
        
        open(file_name, "r") do file           
            while !eof(file)
                line = strip(readline(file))
                state = _define_vtk_state(line)
                iszero(state) && continue
                
                if state == 1
                    try
                        data = split(line, r" +")
                        length(data) == 3 || error("")
                        np = data[2] |> Meta.parse |> eval
                        np isa Int64 && np > 0 || error("")
                        for _ in 1:np
                            line = strip(readline(file))
                            data = split(line, r" +") .|> Meta.parse .|> eval
                            length(data) == 3 || data isa Vector{<:Real} || error("")
                            push!(points, Point2D(data[1], data[2]))
                        end
                    catch
                        error("Ошибка парсинга раздела `POINTS`")
                    end   
                elseif state == 2
                    try
                        data = split(line, r" +")
                        length(data) == 3 || error("")
                        nc = data[2] |> Meta.parse |> eval
                        nc isa Int64 && nc > 0 || error("")
                        for _ in 1:nc
                            line = strip(readline(file))
                            data = split(line, r" +") .|> Meta.parse .|> eval
                            length(data) == 4 || data isa Vector{<:Real} || error("")
                            push!(cells, Cell2D(3, data[2]+1, data[3]+1, data[4]+1))
                        end
                    catch
                        error("Ошибка парсинга раздела `CELLS`")
                    end    
                elseif  state == 3
                    continue
                elseif state == 4
                    try
                        data = split(line, r" +")
                        length(data) == 2 || error("")
                        npd = data[2] |> Meta.parse |> eval
                        npd isa Int64 && npd > 0 || error("")
                    catch
                        error("Ошибка парсинга раздела `POINT_DATA`")
                    end    
                elseif state == 5
                    try
                        line = strip(readline(file))
                        line == "LOOKUP_TABLE default" || error("")
                        for _ in 1:npd
                            line = strip(readline(file))
                            data = line |> Meta.parse |> eval
                            data isa Real || error("")
                            push!(uxs, float(data))
                        end
                    catch
                        error("Ошибка парсинга раздела `SCALARS ux`")
                    end   
                elseif state == 6
                    try
                        line = strip(readline(file))
                        line == "LOOKUP_TABLE default" || error("")
                        for _ in 1:npd
                            line = strip(readline(file))
                            data = line |> Meta.parse |> eval
                            data isa Real || error("")
                            push!(uys, float(data))
                        end
                    catch
                        error("Ошибка парсинга раздела `SCALARS uy`")
                    end    
                end
            end
        end

        return  FiberTaskGridParse(points, cells, uxs, uys)
    end

    function _define_vtk_state(line)::Int64
        contains(line, "POINTS") && return 1
        contains(line, "CELLS") && return 2
        contains(line, "CELL_TYPES") && return 3
        contains(line, "POINT_DATA") && return 4
        contains(line, "SCALARS ux") && return 5
        contains(line, "SCALARS uy") && return 6
        return 0            
    end
end