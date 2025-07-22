module mnemia

# consider the r package JuliaCall

import StatsBase as sb
import LinearAlgebra as la
import Distributions as db
import Plots as pl
import GraphRecipes as gr
import GraphPlot as gp
import ColorTypes as rgba

function transitiveClosure(m::Matrix)
    n = size(m)[1]
    new = copy(m)
    for k in 1:n
        for i in 1:n
            for j in 1:n
                if new[i,k] == 1 && new[k,j] == 1
                    new[i,j] = 1
                end
                if i == j
                    new[i,j] = 1
                end
            end
        end
    end
    new
end

function transitiveReduction(m::Matrix)
    n = size(m)[1]
    new = copy(m)
    for k in 1:n
        for i in 1:n
            for j in 1:n
                if new[i,k] == 1 && new[k,j] == 1 && i != k && j != k
                    new[i,j] = 0
                end
                if i == j
                    new[i,j] = 1
                end
            end
        end
    end
    new
end

function sim(k = 2, s = 5, e = 2, n = 10, p = 0.1, tree = false)
    phis = zeros(s,s,k)
    m = e*s
    thetas = zeros(s,m,k)
    for i in 1:k
        phis[1:s,1:s,i] = sb.sample([0,1],sb.Weights([1-p,p]),s*s)
        for j in 1:s
            for l in 1:s
                if j > l
                    phis[j,l,i] = 0
                end
            end
        end
        if tree
            phis[:,:,i] = makeLikeTree(phis[:,:,i])
        end
        idx = sb.sample(1:s,s,replace=false)
        if s==n
            idx=1:n
        end
        phis[1:s,1:s,i] = phis[idx,idx,i]
        phis[1:s,1:s,i] = transitiveClosure(phis[1:s,1:s,i])
        for j in 1:m
            idx = sb.sample(1:s,1)
            thetas[idx[1],j,i] = 1
        end
        if i == 1
            D = transpose(phis[:,:,i] * thetas[:,:,i])
        else
            D = hcat(D,transpose(phis[:,:,i] * thetas[:,:,i]))
        end
    end
    sidx = sb.sample(1:size(D)[2],n)
    rho = zeros(s,n)
    for i in 1:n
        nidx = mod(sidx[i],s)
        if nidx == 0
            nidx = s
        end
        rho[nidx,i] = 1
    end
    if s==n
        rho=zeros(s,n)
        rho[la.diagind(rho)].=1
        sidx=1:n
    end
    D = D[:,sidx]
    return phis, thetas, rho, D
end

function addNoise(D, mu = 1, sigma = 1)
    norm = db.Normal(0,sigma)
    R = copy(D)
    for i in 1:size(R)[1]
        for j in 1:size(R)[2]
            if D[i,j] == 1
                R[i,j] = mu + db.rand(norm)
            else
                R[i,j] = -mu + db.rand(norm)
            end
        end
    end
    R
end

function score(phi, R, rho)
    S = R * transpose(rho) * transitiveClosure(phi)
    theta = zeros(size(phi)[1],size(R)[1])
    for i in 1:size(R)[1]
        theta[findmax(S[i,:])[2],i] = sb.sample([1],)
    end
    S = theta * S
    Escore=maximum(S;dims=2)
    S = sum(sb.diag(S))
    return S, Escore
end

function addClassNodes(x,y)
    class1=sum(x[:,y.==1],dims=2)
    class0=sum(x[:,y.==0],dims=2)
    class1[class1.>class0].=1
    class0[class0.>class1].=1
    class1[class1.>0].=1
    class0[class0.>0].=1
    xnew=hcat(class1,class0,copy(x))
    xnew
end

function maxSpanTree(R,sparse=false)
    thetamax=zeros(size(R)[1],size(R)[2])
    if sparse
        for i in 1:size(R)[1]
            idx=findmax(R[i,:])
            thetamax[i,idx[2]]=1
        end
    else
        for i in 1:size(R)[2]
            thetamax[findall(R[:,i].>0),i].=1
        end
    end
    tree=zeros(size(R)[2],size(R)[2])
    for i in 1:size(R)[2]
        for j in 1:size(R)[2]
            tree[i,j]=transpose(R[:,i])*thetamax[:,j]
        end
    end
    tree[la.diagind(tree)].=0
    tree[findall(tree.<0)].=0
    phi=zeros(size(R)[2],size(R)[2])
    for i in 1:(size(R)[2]-1)
        if all(tree.<=0)
            break
        end
        maxidx=findmax(tree)
        phi[maxidx[2]]=1
        phi=transitiveClosure(phi)
        tree[maxidx[2]]=0
        phi[la.diagind(phi)].=0
        while maximum(phi+transpose(phi))>1 || sum(phi[:,maxidx[2][2]])>1
            phi[maxidx[2]]=0
            maxidx=findmax(tree)
            phi[maxidx[2]]=1
            phi=transitiveClosure(phi)
            tree[maxidx[2]]=0
            phi[la.diagind(phi)].=0
        end
    end
    phi
end

function makeLikeTree(phi)
    tree=copy(phi)
    treesums=vec(sum(tree,dims=2))
    idx=reverse(sortperm(treesums))
    tree=tree[idx,idx]
    tree[la.tril!(trues(size(tree)),-1)].=0
    tree=tree[sortperm(idx),sortperm(idx)]
    tree[la.diagind(tree)] = repeat([0],size(phi)[1])
    for i in 1:size(tree)[2]
        idx1 = findall(>(0),tree[:,i])
        if size(idx1)[1]>0
            parent = sb.sample(idx1,1)
            tree[:,i] = repeat([0],size(tree)[1])
            tree[parent,i] = [1]
        end
    end
    tree
end

function greedy(R, rho, phi = nothing, tree = false, starts = 10, threads = 1, forbid=nothing, force=nothing)
    if phi === nothing
        phi = transitiveClosure(zeros(size(rho)[1],size(rho)[1]))
    else
        starts=1
    end
    if forbid==nothing
        forbid=copy(phi).*0
    end
    if force==nothing
        force=copy(phi).*0
    end
    bestll=-Inf
    bestphi=copy(phi)
    lls=[]
    phis=[]
    Threads.@threads for run in 1:starts
        if run>1
            phi=transitiveReduction(rand(0:1,size(phi)[1],size(phi)[1]))
            if tree
                phi=makeLikeTree(phi)
            end
            phi[findall(forbid.==1)].=0
            phi[findall(force.==1)].=1
        end
        llold = -Inf
        ll = score(phi,R,rho)[1]
        while ll > llold
            llold = ll
            phiscore = phi*0
            phiscore = convert(Matrix{Float32},phiscore)
            phiscore[:,:] = repeat([-Inf],size(phi)[1]^2)
            for i in 1:size(phi)[1]
                for j in 1:size(phi)[2]
                    if i != j
                        phinew = copy(phi)
                        if tree
                            phinew = transitiveReduction(phinew)
                            if phinew[i,j] == 0
                                phinew[:,j] = repeat([0],size(phinew)[1])
                                phinew[j,i] = 0
                            end
                        end
                        phinew[i,j] = 1 - phinew[i,j]
                        if force[i,j]==1
                            phinew[i,j]=1
                        end
                        if forbid[i,j]==1
                            phinew[i,j]=0
                        end
                        phiscore[i,j] = score(phinew,R,rho)[1]
                    end
                end
            end
            ll = findmax(phiscore)
            if ll[1]>llold
                if tree
                    phi = transitiveReduction(phi)
                    if phi[ll[2]] == 0
                        phi[:,ll[2][2]] = repeat([0],size(phi)[1])
                        phi[ll[2][2],ll[2][1]] = 0
                    end
                end
                phi[ll[2]] = 1 - phi[ll[2]]
                if force[ll[2]]==1
                    phi[ll[2]]=1
                end
                if forbid[ll[2]]==1
                    phi[ll[2]]=0
                end
                phi = transitiveClosure(phi)
            end
            ll = ll[1]
        end
        push!(lls,ll)
        push!(phis,phi)
        if ll>bestll
            bestll=ll
            bestphi=copy(phi)
        end
    end
    return bestphi, bestll, lls, phis
end

function plotNEM(G,tree=false,col=nothing,resolution=(2048,2048))
    n = size(G)[1]
    if col==nothing
        nodecol = [rgba.RGBA{Float64}(1, 0, 0, 0.5)]
        for i in 2:n
            nodecol = hcat(nodecol,rgba.RGBA{Float64}(1, 0, 0, 0.5))
        end
    else
        nodecol=col
    end
    G = transitiveReduction(G)
    G[la.diagind(G)].=0
    # `:spectral`, `:sfdp`, `:circular`, `:shell`, `:stress`, `:spring`, `:tree`, `:buchheim`, `:arcdiagram` or `:chorddiagram`
    method = :stress # circular/shell or stress?
    Gnames = 1:n
    if tree
        G = vcat(transpose(repeat([1],size(G)[2]+1)),hcat(repeat([0],size(G)[1]),G))
        kids = findall(>(1),sum(eachrow(G)))
        G[1,kids] = repeat([0],size(kids)[1])
        G[1,1] = 0
        Gnames = ("root",Gnames...)
        nodecol = hcat(rgba.RGBA{Float64}(1,1,1,1),nodecol)
        method = :buchheim
    end
    nodecol = nodecol[1,:]
    p1 = gr.graphplot(G, names=Gnames, markercolor=nodecol,
    fontsize=10, nodeshape=:circle, method=method)
    pl.plot(p1, layout=(1,1), size=resolution)
end

function resp(phis, R, rho, gamma, complete=false)
    for i in 1:size(phis)[3]
        Rw = R.*transpose(gamma[i,:])
        phi = phis[:,:,i]
        S = Rw * transpose(rho) * phi
        theta = zeros(size(phi)[1],size(R)[1])
        for i in 1:size(R)[1]
            theta[findmax(S[i,:])[2],i] = sb.sample([1],)
        end
        gamma[i,:] = sb.diag(transpose(rho) * phi * theta * R)
    end
    gamma
end

function norm(x)
    xmax = maximum(x)
    maxnum = 2^1023.3
    shrinkfac = log(maxnum)
    x = x .- (xmax - shrinkfac/length(x))
    x
end

function norm2(x)
    x = x/sum(x)
    x
end

function cll(gamma, pi)
    y = mapslices(norm, gamma; dims=1)
    y = exp.(y)
    y = y .* pi
    y = mapslices(norm2, y; dims=1)
    y = y .* (gamma .+ log.(pi))
    ll = sum(mapslices(sum, y; dims=1))
    ll
end

function oll(gamma, pi)
    y = exp.(gamma)
    y = y .* pi
    ll = sum(log.(mapslices(sum,y;dims=1)))
    ll
end

function mnem(R, rho, k = 1, maxiter = 100, tree = false, complete = false)
    s = size(rho)[1]
    phis = sim(k,s,1,1,0.5,tree)[1]
    pi = repeat([1/k],k)
    gamma = ones(k,size(R)[2])/k
    gamma = resp(phis,R,rho,gamma,complete)
    if (complete)
        ll = cll(gamma,pi)
    else
        ll = oll(gamma,pi)
    end
    llold = -Inf
    lls = ll
    iter = 0
    while llold < ll && iter < maxiter
        iter = iter + 1
        llold = ll
        if (complete)
            gamma=mapslices(norm, gamma; dims=1)
        end
        gammaW = exp.(gamma).*pi
        gammaW = gammaW./transpose(transpose(gammaW)*repeat([1],k))
        pi = gammaW*repeat([1],size(gammaW)[2])
        pi = pi./sum(pi)
        for i in 1:k
            Rw = R.*transpose(gammaW[i,:])
            phis[:,:,i] = greedy(Rw, rho, phis[:,:,i], tree)[1]
        end
        gamma = resp(phis,R,rho,gammaW,complete)
        if (complete)
            ll = cll(gamma,pi)
        else
            ll = oll(gamma,pi)
        end
        lls = hcat(lls,ll)
    end
    return phis, lls, gamma, pi, tree, complete
end

function plotMix(x)
    k = size(x[1])[3]
    n = size(x[1])[1]
    gsort = sortperm(x[4],rev=true)
    G = x[1][:,:,gsort[1]]
    nodecol = [rgba.RGBA{Float64}(1, 0, 0, 0.5)]
    for i in 2:n
        nodecol = hcat(nodecol,rgba.RGBA{Float64}(1, 0, 0, 0.5))
    end
    for i in 2:k
        Gtmp = G[:,1:n]*0
        Gtmp2 = G[1:n,:]*0
        G = hcat(G,Gtmp)
        G2 = x[1][:,:,gsort[i]]
        G2 = hcat(Gtmp2,G2)
        G = vcat(G,G2)
        for j in 1:n
            nodecol = hcat(nodecol,rgba.RGBA{Float64}(1-i/k,i/k,i/k, 0.5))
        end
    end
    G = transitiveReduction(G)
    G[la.diagind(G)] = repeat([0],n*k)
    p2 = pl.plot(x[2][1,:], label="model log likelihood",grid=true,gridalpha=0.5,
    ylim=(findmin(x[2][1,:])[1]-1,findmax(x[2][1,:])[1]*1.1),
    plot_title="convergence",plot_titlefontsize=10,legend=:bottomright)
    gamma = x[3]
    gammaW = exp.(gamma).*x[4]
    gammaW = gammaW./transpose(transpose(gammaW)*repeat([1],k))
    p3 = pl.histogram(gammaW[gsort[1],:],fillcolor=rgba.RGBA{Float64}(1, 0, 0, 0.5),
    labels=transpose(round(x[4][gsort[1]],digits=2)),bins=20,
    plot_title="sample posterior",plot_titlefontsize=10,legend=:top)
    for i in 2:k
        pl.histogram!(gammaW[gsort[i],:],fillcolor=rgba.RGBA{Float64}(1-i/k,i/k,i/k, 0.5),
        labels=transpose(round(x[4][gsort[i]],digits=2)),bins=20)
    end
    method = :stress # circular/shell or stress?
    if k > 1
        s = floor(Int,k/2)
        G2 = G[(s*n+1):(n*k),(s*n+1):(n*k)]
        G = G[1:(s*n),1:(s*n)]
        Gnames = repeat(1:n,s)
        G2names = repeat(1:n,(k-s))
        nodecol2 = nodecol[:,(s*n+1):(k*n)]
        nodecol = nodecol[:,1:(s*n)]
        if x[5]
            G = vcat(transpose(repeat([1],size(G)[2]+1)),hcat(repeat([0],size(G)[1]),G))
            kids = findall(>(1),sum(eachrow(G)))
            G[1,kids] = repeat([0],size(kids)[1])
            G[1,1] = 0
            G2 = vcat(transpose(repeat([1],size(G2)[2]+1)),hcat(repeat([0],size(G2)[1]),G2))
            kids = findall(>(1),sum(eachrow(G2)))
            G2[1,kids] = repeat([0],size(kids)[1])
            G2[1,1] = 0
            Gnames = ("root",Gnames...)
            G2names = ("root",G2names...)
            nodecol = hcat(rgba.RGBA{Float64}(1,1,1,1),nodecol)
            nodecol2 = hcat(rgba.RGBA{Float64}(1,1,1,1),nodecol2)
            method = :buchheim
        end
        nodecol = nodecol[1,:]
        nodecol2 = nodecol2[1,:]
        p0 = gr.graphplot(G, names=Gnames, markercolor=nodecol,
        fontsize=10, nodeshape=:circle, method=method)
        p1 = gr.graphplot(G2, names=G2names, markercolor=nodecol2,
        fontsize=10, nodeshape=:circle, method=method)
        pl.plot(p0, p1, p2, p3, layout=(2,2))
    else
        Gnames = repeat(1:n,k)
        if x[5]
            G = vcat(transpose(repeat([1],size(G)[2]+1)),hcat(repeat([0],size(G)[1]),G))
            kids = findall(>(1),sum(eachrow(G)))
            G[1,kids] = repeat([0],size(kids)[1])
            G[1,1] = 0
            Gnames = ("root",Gnames...)
            nodecol = hcat(rgba.RGBA{Float64}(1,1,1,1),nodecol)
            method = :buchheim
        end
        nodecol = nodecol[1,:]
        p1 = gr.graphplot(G, names=repeat(1:n,k), markercolor=nodecol,
        fontsize=10, nodeshape=:circle, method=method)
        pl.plot(p1, p2, p3, layout=(1,3))
    end
end

function transformData(D, fpr = 0.01, fnr = 0.1)
    R = D
    idx1 = findall(==(1),D)
    R[idx1] = repeat([log((1-fnr)/fpr)],size(idx1)[1])
    idx0 = findall(==(0),D)
    R[idx0] = repeat([log(fnr/(1-fpr))],size(idx0)[1])
    R
end

# https://turinglang.org/docs/getting-started/

function samplePhis(phis)
    phisnew=copy(phis)
    for i in 1:length(phis)
        move=rand((1:2))
        if move==1
            ycord=rand((1:size(phis[i])[1]))
            xcord=rand((1:size(phis[i])[2]))
            while xcord==ycord
                xcord=rand((1:size(phis[i])[2]))
            end
            phisnew[1]=copy(phis[i])
            phisnew[1][ycord,xcord]=1-phisnew[1][ycord,xcord]
        end
        if move==2
            ycord=rand((1:size(phis[i])[1]))
            xcord=rand((1:size(phis[i])[2]))
            while xcord==ycord
                xcord=rand((1:size(phis[i])[2]))
            end
            phisnew[1]=copy(phis[i])
            phisnew[1][ycord,xcord]=0
            phisnew[1][xcord,ycord]=1
        end
    end
    phisnew
end

function sampleThetas(phis,R,rho,gamma)
    thetas=[]
    for i in 1:size(phis)[3]
        Rw = R.*transpose(gamma[i,:])
        phi = phis[:,:,i]
        S = Rw * transpose(rho) * phi
        theta = zeros(size(phi)[1],size(R)[1])
        for i in 1:size(R)[1]
            theta[findmax(S[i,:])[2],i] = sb.sample([1],)
        end
        gamma[i,:] = sb.diag(transpose(rho) * phi * theta * R)
        push!(thetas,theta)
    end
    thetas
end

function mcmc(R, rho, k = 1, iters = 1000, burnin = 100, tree = false, complete = false)
    
end



end