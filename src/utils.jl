module Utils
    export _define_conforming_point,
        _define_period

    function _define_conforming_point(gminx::Float64, gmaxx::Float64, gminy::Float64, gmaxy::Float64, x::Float64, y::Float64)::Tuple{Float64, Float64, Bool, Bool}
        L = gmaxx - gminx
        H = gmaxy - gminy
        while x < gminx || y < gminy
            x += 100 * L
            y += 100 * H
        end
        flagx = iseven(div(x - gminx, L))
        flagy = iseven(div(y - gminy, H))
        newx = flagx ? mod(x - gminx, L) + gminx : L - mod(x - gminx, L) + gminx
        newy = flagy ? mod(y - gminy, H) + gminy : H - mod(y - gminy, H) + gminy

        return newx, newy, flagx, flagy
    end

    function _define_period(x, x_min, L)
        x_period = if x < x_min + L && x >= x_min
            1
        elseif x >= x_min + L
            x_period = 1
            while x >= x_min + L
                x -= L 
                x_period += 1
            end
            x_period
        elseif x < x_min
            x_period = 0
            while x < x_min
                x += L
                x_period -= 1
            end
            x_period
        end
        return x_period
    end

    _define_period(x, y, x_min, L, y_min, H) = (_define_period(x, x_min, L), _define_period(y, y_min, H))
end