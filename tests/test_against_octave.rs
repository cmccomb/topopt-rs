use mocktave::eval;
use nalgebra::DMatrix;

fn mocktave_top88(
    nelx: usize,
    nely: usize,
    volfrac: f64,
    penalMax: f64,
    rmin: f64,
) -> DMatrix<f64> {
    // This code adapted from here: https://github.com/blademwang11/Topopt/blob/master/top88.m
    let script = format!(
        "
    function xPhys = top88(nelx,nely,volfrac,penalMax,rmin)
        E0 = 1;
        Emin = 1e-9;
        nu = 0.3;
        penal = 0.96;
        A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
        A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
        B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
        B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
        KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
        nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
        edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
        edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
        iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
        jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
        F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
        U = zeros(2*(nely+1)*(nelx+1),1);
        fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
        alldofs = [1:2*(nely+1)*(nelx+1)];
        freedofs = setdiff(alldofs,fixeddofs);
        iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
        jH = ones(size(iH));
        sH = zeros(size(iH));
        k = 0;
        for i1 = 1:nelx
          for j1 = 1:nely
            e1 = (i1-1)*nely+j1;
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
              for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
              endfor
            endfor
          endfor
        endfor
        H = sparse(iH,jH,sH);
        Hs = sum(H,2);
        x = repmat(volfrac,nely,nelx);
        xPhys = x;
        loop = 0;
        change = 1;
        while (change > 0.01)
          loop = loop + 1;
          penal = min(penalMax, penal + 0.04);
          %% FE-ANALYSIS
          sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
          K = sparse(iK,jK,sK); K = (K+K')/2;
          tic; U(freedofs) = K(freedofs,freedofs)\\F(freedofs); toc;
          ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
          c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
          dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
          dv = ones(nely,nelx);
          dc(:) = H*(dc(:)./Hs);
          dv(:) = H*(dv(:)./Hs);
          l1 = 0; l2 = 1e9; move = 0.2;
          while ((l2-l1)/(l1+l2) > 1e-3)
            lmid = 0.5*(l2+l1);
            xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
            xPhys(:) = (H*xnew(:))./Hs;
            if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; endif
          endwhile
          change = max(abs(xnew(:)-x(:)));
          x = xnew;
        endwhile
    endfunction

    z = top88({nelx},{nely},{volfrac},{penalMax},{rmin});

    "
    );

    let y = mocktave::eval(&script).get_matrix_named("z").expect("Wtf?");

    DMatrix::from_fn(nelx, nely, |i, j| y[i][j])
}

#[cfg(test)]
mod tests {
    use core::assert_eq;

    #[test]
    fn silly_simple() {
        let settings = topopt::Settings::new(30, 10, 0.5)
            .with_left_bc(true, false)
            .with_bottom_right_bc(false, true)
            .with_top_left_load(0.0, -1.0)
            .with_penalty_weight(3.0)
            .with_filter_radius(1.5);

        let x = topopt::solve(settings);
        // let y = super::mocktave_top88(30, 10, 0.5, 3.0, 1.5);
        // assert_eq!(x, y)
    }
}
