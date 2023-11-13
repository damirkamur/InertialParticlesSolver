@info "Import IntertialParticlesSolver"
using IntertialParticlesSolver
@info "Start main"

function main()
    Sts = [10^i for i=LinRange(-3, 2, 200)]
    N = 10000

    task_name = "periodic_geom_f_53"
    max_period = 10
    @info "Загрузка задачи"
    p = Parser(task_name)
    particles = [Particle(0.0, i, 0.1, 0.0) for i ∈ LinRange(0, 1, N)]
    @info "Генерация задачи"
    task = task_from_parse(
        p, 
        deepcopy(particles); 
        tend = 100, 
        max_period = max_period,
        St = 0.02, 
        N = 1000, 
        M = 1000, 
        delta = 0.0, 
        periodic = true
    )
    for p in particles
        (vx0, vy0) = get_value(task.grid, p.x0, p.y0)
        p.vx0 = vx0
        p.vy0 = vy0
    end
    if !isdir("data")
        mkdir("data")
    end
    folder = "data/" * task_name * " N = $N"
    if !isdir(folder)
        mkdir(folder)
    end
    task_name = folder * "/data"

    open("$task_name.txt", "w") do file 
        println(file, "# N")
        println(file, "# NSt's")
        println(file, "# St's")
        println(file, "# max_period")
        println(file, "# (ESt_1) term_1 term_2 ... term_max_period")
        println(file, "# ...")
        println(file, "# (ESt_NSt's) term_1 term_2 ... term_max_period")
        
        println(file, "$N")
        println(file, "$(length(Sts))")
        println(file, join(Sts, " "))
        println(file, "$max_period")
    end
    @info "Начало цикла"

    for (index, St) in enumerate(Sts)
        @info "Задача $index) St = $St"
        task.particles = [Particle(0.0, i, 0.1, 0.0) for i ∈ LinRange(0, 1, N)]
        for p in task.particles
            (vx0, vy0) = get_value(task.grid, p.x0, p.y0)
            p.vx0 = vx0
            p.vy0 = vy0
        end
        task.St = St
        task.terminated = 0
        
        @info "Решение"
        t1 = time()
        solve(task)
        t2 = time()
        @info "Решение заняло $(round(t2 - t1; digits = 2)) сек."
        @info "Сбор информации"
        ti = TaskInfo(task)
        @info "Запись информации в файл"
        open("$task_name.txt", "a") do file
            for period = 1:max_period
                print(file, "$(get_xperiod_catch(ti, period)) ")
            end
            println(file)
        end
        GC.gc()

        @info "Рисование картинок"
        draw(task, 3000; 
            file_name="$task_name" * " $index) St = $(ti.St) N = $(ti.N) terminated = $(ti.terminated).png", 
            dy_left = 0, 
            dy_right = 0, 
            dx_left = 0, 
            dx_right = 2, 
            colorize_terminated=true, 
            dark_theme=true
        )
        # draw(task, 3000; 
        #     file_name="$task_name" * " $index) St = $(ti.St) N = $(ti.N) terminated = $(ti.terminated).svg", 
        #     dy_left = 0, 
        #     dy_right = 0, 
        #     dx_left = 0, 
        #     dx_right = 2, 
        #     colorize_terminated=true, 
        #     dark_theme=true
        # )
        # draw(task, 3000; 
        #     file_name="$task_name" * " $index) St = $(ti.St) N = $(ti.N) terminated = $(ti.terminated).eps", 
        #     dy_left = 0, 
        #     dy_right = 0, 
        #     dx_left = 0, 
        #     dx_right = 2, 
        #     colorize_terminated=true, 
        #     dark_theme=true
        # )

    end
end

main()

