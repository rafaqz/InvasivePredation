simplify!(ax; kw...) = (hidedecorations!(ax; kw...); hidespines!(ax))

function label!(fig, text, i, j;
    fontsize=20
)
    text_ax = Axis(fig[i, j])
    xlims!(text_ax, (0, 1))
    ylims!(text_ax, (0, 1))
    simplify!(text_ax)
    text!(text_ax, 0.0, 0.0; text=string(text))

    return text_ax
end
