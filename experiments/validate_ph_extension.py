#!/usr/bin/env python3
"""Discriminating nonlinear/open PH validation of the geometric time extension.

Uses the unmodified 2012 mixed spatial construction in the v12 verifier,
expressed in orthonormal Legendre coordinates. Run from any directory.
Outputs are placed in experiments/results and experiments/figures.
"""
from __future__ import annotations

import os
for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                  "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_variable, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/ph-extension-mpl")

import csv
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from numpy.polynomial.legendre import leggauss, legvander

from verify_transmission_line_v12 import shifted_legendre_line, structural_residuals
from ph_gauges import PHModel, scheme, common_gauge

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"


def quartic_line(N=6, alpha=2.0, damping=0.0, forcing=None):
    """q_t=-d_z p, p_t=-d_z(q+alpha*q^3), plus modal PH losses/ports.

    Unit length, capacitance and inductance make the zero-state energy
    metric Q0=I in the orthonormal modal chart. Nonlinear constitutive
    efforts are projected gradients of the integrated Hamiltonian.
    """
    data = shifted_legendre_line(N,1.,1.,1.)
    xi,w = leggauss(max(2*N+2,32))
    z,w = (xi+1)/2,w/2
    V = legvander(xi,N-1)*np.sqrt(2*np.arange(N)+1)
    def H(x):
        q,p = V@x[:N],V@x[N:]
        return float(w@(.5*(q*q+p*p)+.25*alpha*q**4))
    def grad(x):
        q,p = V@x[:N],V@x[N:]
        return np.r_[V.T@(w*(q+alpha*q**3)),V.T@(w*p)]
    def hess(x):
        q = V@x[:N]
        out = np.zeros((2*N,2*N))
        out[:N,:N] = V.T@((w*(1+3*alpha*q*q))[:,None]*V)
        out[N:,N:] = V.T@(w[:,None]*V)
        return out
    model = PHModel(H,grad,hess,data.interconnection,damping*np.eye(2*N),
                    data.input_map,data.feedthrough,forcing)
    model.N,model.spatial,model.alpha = N,data,alpha
    model.z,model.w,model.V = z,w,V
    return model


def initial(model,amplitude=.35):
    return np.r_[model.V.T@(model.w*amplitude*np.sin(np.pi*model.z/2)),
                 np.zeros(model.N)]


def spatial_theta(model):
    """Noncentral alternating form derived from the spatial differential.

    D M is the derivative restricted to active Legendre modes. The form
    below is not chosen by random entries or by fitting a reference flow.
    Its Frobenius normalization is intrinsic for this Q0=I modal chart.
    """
    active = model.spatial.derivative @ model.spatial.mixed_mass
    zero = np.zeros_like(active)
    raw = np.block([[zero,active],[-active.T,zero]])
    return raw/np.linalg.norm(raw,"fro")


def vector_field(model,t,x):
    u = np.zeros(2) if model.u is None else model.u(t)
    return (model.J-model.R)@model.grad(x)+model.B@u


def reference_pair(model,x0,T,times=None):
    fun = lambda t,x: vector_field(model,t,x)
    tight = solve_ivp(fun,(0,T),x0,method="DOP853",rtol=2.5e-13,atol=2.5e-15,
                      t_eval=times,dense_output=True)
    loose = solve_ivp(fun,(0,T),x0,method="DOP853",rtol=1e-11,atol=1e-13)
    if not tight.success or not loose.success:
        raise RuntimeError("Reference solve failed")
    discrepancy = float(np.linalg.norm(tight.y[:,-1]-loose.y[:,-1]))
    return tight,discrepancy


def integrate(model,x0,T,h,s,q,gauged=False):
    n = int(round(T/h))
    if abs(n*h-T)>1e-12:
        raise ValueError("T must be an integer multiple of h")
    integrator = scheme(s=s,q=q,gauge=common_gauge(spatial_theta(model)) if gauged else None)
    x = x0.copy()
    states,energy,diags = [x.copy()],[model.H(x)],[]
    for i in range(n):
        x,d = integrator.step(model,x,i*h,h)
        states.append(x.copy());energy.append(model.H(x));diags.append(d)
    ts = np.linspace(0,T,n+1)
    states,energy = np.asarray(states),np.asarray(energy)
    work = np.r_[0.,np.cumsum([d["work"] for d in diags])]
    diss = np.r_[0.,np.cumsum([d["diss"] for d in diags])]
    record = {"N":model.N,"T":T,"h":h,"s":s,"q":q,"gauged":gauged,"steps":n,
              "energy_initial":float(energy[0]),"energy_final":float(energy[-1]),
              "max_relative_energy_change":float(np.max(np.abs(energy-energy[0]))/max(abs(energy[0]),1e-30)),
              "max_cumulative_power_defect":float(np.max(np.abs(energy-energy[0]-work+diss))),
              "total_work":float(work[-1]),"total_dissipation":float(diss[-1]),
              "max_positive_energy_increment":float(max(0.,np.max(np.diff(energy)))),
              "max_iterations":int(max(d["iterations"] for d in diags)),
              "min_step_dissipation":float(min(d["diss"] for d in diags))}
    for key in ("balance","chain_defect","gauge_power","structure_defect",
                "residual_power_certificate","certificate_closure","stage_residual"):
        record["max_abs_"+key] = float(max(abs(d[key]) for d in diags))
    record["max_effort_correction"] = float(max(np.linalg.norm(d["e"]-d["G"]) for d in diags))
    return record,ts,states,energy,work,diss,diags


def structural_checks():
    records=[]
    for N in (4,6,8):
        model=quartic_line(N)
        x=initial(model)
        theta=spatial_theta(model)
        A=model.J
        direction=np.linspace(-1,1,2*N)
        eps=1e-6
        gd=abs((model.H(x+eps*direction)-model.H(x-eps*direction))/(2*eps)
                -model.grad(x)@direction)
        hd=np.linalg.norm((model.grad(x+eps*direction)-model.grad(x-eps*direction))/(2*eps)
                          -model.hess(x)@direction,np.inf)
        row={"N":N,**{k:float(v) for k,v in structural_residuals(model.spatial).items()},
             "theta_skew_defect":float(np.linalg.norm(theta+theta.T)),
             "theta_metric_norm":float(np.linalg.norm(theta,"fro")),
             "generator_gauge_noncommutation":float(np.linalg.norm(model.J@theta@A-A@model.J@theta)),
             "directional_gradient_error":float(gd),"directional_hessian_error":float(hd)}
        assert row["stokes"]<2e-10 and row["directional_gradient_error"]<1e-8
        assert row["directional_hessian_error"]<1e-7 and row["generator_gauge_noncommutation"]>1
        records.append(row)
    return records


def temporal_order():
    model=quartic_line(6);x0=initial(model);T=1.
    ref,refgap=reference_pair(model,x0,T)
    rows=[]
    for s,q in ((1,2),(2,4)):
        for gauged in (False,True):
            previous=None
            for h in (.1,.05,.025,.0125,.00625):
                r,_,xx,_,_,_,_=integrate(model,x0,T,h,s,q,gauged)
                err=float(np.linalg.norm(xx[-1]-ref.y[:,-1]))
                r.update(name=f"PG{2*s}"+(" + spatial gauge" if gauged else ""),
                         final_L2_error=err,observed_order=None if previous is None else float(np.log2(previous/err)))
                assert r["max_relative_energy_change"]<1e-10
                assert r["max_abs_certificate_closure"]<1e-12
                rows.append(r);previous=err
            assert rows[-1]["observed_order"]>2*s-.3
    assert refgap<1e-10
    return {"reference_tolerance_discrepancy_L2":refgap,"rows":rows}


def open_damped():
    forcing=lambda t:np.array([.06*np.sin(1.7*t),.04*np.sin(2.3*t)])
    model=quartic_line(8,damping=.12,forcing=forcing)
    x0=initial(model,.12);T,h=4.,.04
    times=np.linspace(0,T,int(T/h)+1)
    ref,refgap=reference_pair(model,x0,T,times)
    rows=[];plotdata=None;previous_error={False:None,True:None}
    for h,gauged in ((step,flag) for step in (.04,.02,.01,.005,.0025) for flag in (False,True)):
        r,ts,xx,en,work,diss,dd=integrate(model,x0,T,h,2,4,gauged)
        r["final_L2_error"]=float(np.linalg.norm(xx[-1]-ref.y[:,-1]))
        old=previous_error[gauged]
        r["observed_order"]=None if old is None else float(np.log2(old/r["final_L2_error"]))
        previous_error[gauged]=r["final_L2_error"]
        # Recover physical boundary efforts from the original rectangular data.
        N=model.N;data=model.spatial
        left=np.vstack([data.mixed_mass.T,data.left_trace])
        right=np.vstack([data.mixed_mass.T,data.right_trace])
        physical_output_defect=0.;density_balance_defect=0.
        for d in dd:
            for e,u,y,F in zip(d["e"],d["u"],d["y"],d["F"]):
                a=np.linalg.solve(left,np.r_[e[:N],u[0]])
                b=np.linalg.solve(right,np.r_[e[N:],u[1]])
                yp=np.array([data.left_trace@b,-data.right_trace@a])
                physical_output_defect=max(physical_output_defect,np.max(np.abs(yp-y)))
                # On L=1 the zeroth orthonormal modal coefficient is the integral.
                loss=model.R@e
                density_balance_defect=max(density_balance_defect,
                    abs(F[0]-(yp[0]-u[1])+loss[0]),
                    abs(F[N]-(u[0]+yp[1])+loss[N]))
        r["physical_boundary_output_max_defect"]=float(physical_output_defect)
        r["integrated_density_balance_max_defect"]=float(density_balance_defect)
        assert r["max_cumulative_power_defect"]<2e-10 and r["min_step_dissipation"]>=-1e-15
        assert physical_output_defect<1e-11 and density_balance_defect<1e-10
        rows.append(r)
        print(f"Forced PG4 gauged={gauged} h={h:g}: error={r['final_L2_error']:.6e}, order={r['observed_order']}",flush=True)
        if gauged:
            plotdata=(ts,xx,en,work,diss,ref.sol(ts).T)
    for gauged in (False,True):
        final_row=[r for r in rows if r["gauged"]==gauged][-1]
        assert final_row["observed_order"]>3.7, "Forced nonlinear refinement must resolve fourth order"
    unforced=quartic_line(6,damping=.12)
    rd,*_=integrate(unforced,initial(unforced),3.,.05,2,4,True)
    assert rd["max_positive_energy_increment"]<1e-12
    return {"reference_tolerance_discrepancy_L2":refgap,"forced_rows":rows,
            "unforced_gauged":rd},plotdata


def insufficient_quadrature():
    model=quartic_line(6,alpha=3.);x0=initial(model,.6)
    records=[]
    for q in (1,2):
        r,*_=integrate(model,x0,2.,.1,1,q,False)
        records.append(r)
    assert records[0]["max_abs_chain_defect"]>1000*max(records[1]["max_abs_chain_defect"],1e-18)
    assert records[0]["max_relative_energy_change"]>1e-7
    return {"description":"Same s=1 path: q=1 midpoint does not integrate the quartic chain rule; q=2 does.",
            "rows":records}


def linear_discrimination():
    model=quartic_line(6,alpha=0.)
    dim=2*model.N;Q=model.hess(np.zeros(dim));A=model.J@Q
    omega=np.linalg.inv(model.J)
    records=[]
    for s,q in ((1,1),(2,2)):
        for gauged in (False,True):
            integrator=scheme(s,q,gauge=common_gauge(spatial_theta(model)) if gauged else None)
            columns=[]
            for j in range(dim):
                value,_=integrator.step(model,np.eye(dim)[:,j],0.,.1)
                columns.append(value)
            M=np.column_stack(columns)
            record={"s":s,"q":q,"gauged":gauged,"h":.1,
                    "energy_orthogonality_defect":float(np.linalg.norm(M.T@Q@M-Q,"fro")/np.linalg.norm(Q,"fro")),
                    "symplectic_defect":float(np.linalg.norm(M.T@omega@M-omega,"fro")/np.linalg.norm(omega,"fro")),
                    "normalized_generator_commutator":float(np.linalg.norm(M@A-A@M,"fro")/(np.linalg.norm(M,"fro")*np.linalg.norm(A,"fro"))),
                    "exact_flow_relative_error":float(np.linalg.norm(M-expm(.1*A),"fro")/np.linalg.norm(expm(.1*A),"fro"))}
            assert record["energy_orthogonality_defect"]<1e-11
            if gauged:
                assert record["normalized_generator_commutator"]>1e-9
            else:
                assert record["normalized_generator_commutator"]<1e-12
            records.append(record)
    # Energy alone cannot detect unresolved phase. Use an exact semidiscrete mode.
    vals,vecs=np.linalg.eig(A)
    ix=np.where(vals.imag>1e-8)[0]
    selected=ix[np.argsort(vals.imag[ix])[3]]
    frequency=float(vals[selected].imag)
    a,b=vecs[:,selected].real,vecs[:,selected].imag
    normalization=np.sqrt(a@Q@a);a,b=a/normalization,b/normalization
    T,h=4.,.2
    r,ts,xx,en,*_=integrate(model,a,T,h,2,2,False)
    cosine=xx@(Q@a)/(a@Q@a);sine=-xx@(Q@b)/(b@Q@b)
    phase=np.unwrap(np.arctan2(sine,cosine))
    exact=np.cos(frequency*T)*a-np.sin(frequency*T)*b
    r.update(frequency=frequency,omega_h=frequency*h,
             phase_lag_rad=float(frequency*T-phase[-1]),
             final_L2_error=float(np.linalg.norm(xx[-1]-exact)))
    assert r["max_relative_energy_change"]<1e-10 and abs(r["phase_lag_rad"])>.1
    return {"one_step_maps":records,"energy_without_trajectory_accuracy":r}


def make_figures(order,plotdata):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "axes.spines.top":False,"axes.spines.right":False})
    FIGURES.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(6.5,4))
    for name in dict.fromkeys(r["name"] for r in order["rows"]):
        rows=[r for r in order["rows"] if r["name"]==name]
        ax.loglog([r["h"] for r in rows],[r["final_L2_error"] for r in rows],"o-",label=name)
    ax.set(xlabel="Time step",ylabel="L2 error against tight DOP853",title="Nonlinear quartic wave: temporal order")
    ax.legend();ax.grid(alpha=.25);fig.tight_layout()
    fig.savefig(FIGURES/"ph_extension_order.png",dpi=170);plt.close(fig)
    ts,xx,en,work,diss,reference=plotdata
    fig,axes=plt.subplots(1,3,figsize=(12,3.7))
    axes[0].plot(ts,en,label="Stored energy");axes[0].plot(ts,work,label="Input work");axes[0].plot(ts,diss,label="Dissipation")
    axes[0].set(xlabel="Time",title="Forced dissipative line");axes[0].legend()
    axes[1].plot(ts,en-en[0]-work+diss);axes[1].set(xlabel="Time",title="Cumulative power defect")
    axes[1].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    axes[2].plot(ts,xx[:,0],label="Gauged PG4");axes[2].plot(ts,reference[:,0],"--",label="DOP853")
    axes[2].set(xlabel="Time",ylabel="Integral of q",title="Independent trajectory check");axes[2].legend()
    for ax in axes:ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(FIGURES/"ph_extension_open_damped.png",dpi=170);plt.close(fig)


def execute():
    RESULTS.mkdir(exist_ok=True)
    data={"environment":{"python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,
                         "platform":platform.platform(),"float":"float64"},
          "scope":"Nonlinear finite-dimensional PH systems obtained from the 2012 spatial reduction. Temporal order is at fixed N; this is not a simultaneous PDE convergence theorem.",
          "gauge":"K=h^(2s)(b b^T) tensor Theta; Theta=[[0,D M],[-(D M)^T,0]] normalized in the Q0=I metric.",
          "constitutive_quadrature":"max(2N+2,32) spatial Gauss points; exact for the quartic energy density in exact arithmetic."}
    data["spatial_and_constitutive_checks"]=structural_checks()
    print("Spatial and constitutive checks passed",flush=True)
    data["nonlinear_temporal_order"]=temporal_order()
    print("Nonlinear temporal orders passed",flush=True)
    data["open_dissipative"],plotdata=open_damped()
    print("Open/dissipative checks passed",flush=True)
    data["insufficient_quadrature"]=insufficient_quadrature()
    data["linear_discriminating_checks"]=linear_discrimination()
    data["all_assertions_passed"]=True
    (RESULTS/"ph_extension_validation.json").write_text(json.dumps(data,indent=2,allow_nan=False),encoding="utf-8")
    rows=data["nonlinear_temporal_order"]["rows"]+data["open_dissipative"]["forced_rows"]
    fields=sorted(set().union(*(row.keys() for row in rows)))
    with (RESULTS/"ph_extension_validation.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    make_figures(data["nonlinear_temporal_order"],plotdata)
    print(json.dumps({"all_assertions_passed":True,
          "final_orders":[r for r in data["nonlinear_temporal_order"]["rows"] if r["h"]==.00625],
          "forced":data["open_dissipative"]["forced_rows"],
          "linear":data["linear_discriminating_checks"]},indent=2))


if __name__=="__main__":
    execute()
