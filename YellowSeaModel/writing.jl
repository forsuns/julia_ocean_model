	function WriteRho(nx,ny,nz,rho)
			fname = "./output/rho.csv"
			csv1 = open(fname,"a+")
			@inbounds for j = 2:ny+1, i = 2:nx+1
					write(csv1, [string.(i),",",string.(j),",",string.(rho[i,j,2:nz+1],","),"\n"])
			end
			close(csv1)
	end
	
	function WriteU(nx,ny,nz,u,umm)
			fname = "./output/u.csv"
			csv2 = open(fname,"a+")
			@inbounds for j = 2:ny+1, i = 2:nx+1
				umm[2:nz+1] = 0.5*(u[i,j,2:nz+1] + u[i-1,j,2:nz+1])
				write(csv2, [string.(i),",",string.(j),",",string.(umm[2:nz+1],","),"\n"])
			end
			close(csv2)
	end
	
	function WriteV(nx,ny,nz,v,vmm)
			fname = "./output/v.csv"
			csv3 = open(fname,"a+")
			@inbounds for j = 2:ny+1, i = 2:nx+1
			vmm[2:nz+1] = 0.5*(v[i,j,2:nz+1] + v[i,j-1,2:nz+1])
			write(csv3, [string.(i),",",string.(j),",",string.(vmm[2:nz+1],","),"\n"])	
			end
			close(csv3)
	end
	
	function WriteW(nx,ny,nz,w)
			fname = "./output/w.csv"
			csv4 = open(fname,"a+")
			@inbounds for j = 2:ny+1, i = 2:nx+1
			write(csv4, [string.(i),",",string.(j),",",string.(w[i,j,2:nz+1],","),"\n"])
			end
			close(csv4)
	end
	
	function WriteT(nx,ny,nz,t)
			fname = "./output/t.csv"
			csv5 = open(fname,"a+")
			@inbounds for j = 2:ny+1, i = 2:nx+1
			write(csv5, [string.(i),",",string.(j),",",string.(t[i,j,2:nz+1],","),"\n"])
			end
			close(csv5)
	end
	
	function WriteS(nx,ny,nz,s)
			fname = "./output/s.csv"
			csv6 = open(fname,"a+")
			@inbounds for j = 2:ny+1, i = 2:nx+1
			write(csv6, [string.(i),",",string.(j),",",string.(s[i,j,2:nz+1],","),"\n"])
			end
			close(csv6)
	end
	
	function WriteA(nx,ny,nz,age)
			fname = "./output/age.csv"
			csv7 = open(fname,"a+")
			@inbounds for j = 2:ny+1, i = 2:nx+1
			write(csv7, [string.(i),",",string.(j),",",string.(age[i,j,2:nz+1],","),"\n"])
			end
			close(csv7)
	end
	
	function WriteQ(nx,ny,nz,q)
			fname = "./output/q.csv"
			csv8 = open(fname,"a+")
			for j = 2:ny+1, i = 2:nx+1
					write(csv8, [string.(i),",",string.(j),",",string.(q[i,j,2:nz+1],","),"\n"])
			end
			close(csv8)
	end