module triangle

using LaTeXStrings
using PyPlot

function bin(x, grid)
  i = 0
  while true
    i += 1
    if x <= grid[i]+step(grid)/2 || i >= length(grid)
      break
    end
  end
  return i
end

tiime = range(0, 1000, step=1e-4)
space = range(-1,1, length=51)

dX = sqrt(step(tiime)) .* randn(length(tiime))
dY = sqrt(step(tiime)) .* randn(length(tiime))

X, Y = zeros(length(tiime)), zeros(length(tiime))

y_fixed = bin(0.0, space)
density_x = zeros(length(space))
for i in 2:length(tiime)
    X[i], Y[i] = X[i-1] + dX[i], Y[i-1] + dY[i]
    if X[i] < -1.0 || Y[i] < -X[i] || X[i] > 1.0
        X[i] = X[i-1]
    end
    if Y[i] < -1.0 || Y[i] < -X[i] || Y[i] > 1.0
        Y[i] = Y[i-1]
    end
    # if bin(Y[i], space) == y_fixed
        bin_x = bin(X[i], space)
        dx = bin_x == 1 || bin_x == length(space) ? step(space) / 2 : step(space)
        density_x[bin_x] += step(tiime)/dx / tiime[end]
    # end
end

if true
    fig, ax = plt.subplots(1,2, figsize=(6,3))
    # plt.scatter(X, Y, s=0.1, c=range(-1,1,length=length(tiime)))
    ax[1].plot(X[1:80001], Y[1:80001], lw=0.3, c="black", alpha=0.4)
    ax[1].set_xlim(-1,1)
    ax[1].set_ylim(-1,1)
    # ax[1].axhline(0, color="black", ls="--", lw=0.3)
    ax[1].plot(range(-1,1,length=10), -range(-1,1,length=10), color="black", lw=0.8)
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[1].set_xlabel("Latent DOF")
    ax[1].set_ylabel("Projected DOF")

    ax[2].plot(space, density_x)
    ax[2].set_ylim(0)
    ax[2].set_xlim(-1,1)
    ax[2].set_yticks([])

    ax[2].plot(X[1:100001], tiime[1:100001], lw=0.8, c="black")
    # ax[2].set_xlim(-1,1)
    # ax[2].set_ylim(0)
    # ax[2].set_yticks([])

    fig.subplots_adjust(wspace=-0.00, hspace=0.0)
    plt.savefig("2d-triangular.pdf")
    plt.close()
end

end
