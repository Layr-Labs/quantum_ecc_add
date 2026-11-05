//! Arithmetic that keeps the provided offset in classical bit registers.
use super::*;

fn bit_op(b: &mut B, kind: OperationType, target: BitId, cond: Option<BitId>) {
    let mut op=Op::empty();op.kind=kind;op.c_target=target;
    if let Some(c)=cond {op.c_condition=c;}
    b.ops.push(op);
}
fn constant(b: &mut B, value: bool) -> BitId {
    let out=b.alloc_bit();
    bit_op(b,if value {OperationType::BitStore1} else {OperationType::BitStore0},out,None);
    out
}
fn copy(b: &mut B, a: BitId) -> BitId {
    let out=constant(b,false);bit_op(b,OperationType::BitStore1,out,Some(a));out
}
fn not(b: &mut B, a: BitId) -> BitId {
    let out=constant(b,true);bit_op(b,OperationType::BitInvert,out,Some(a));out
}
fn xor(b: &mut B, a: BitId, c: BitId) -> BitId {
    let out=copy(b,a);bit_op(b,OperationType::BitInvert,out,Some(c));out
}
fn and(b: &mut B, a: BitId, c: BitId) -> BitId {
    let out=copy(b,a);let nc=not(b,c);
    bit_op(b,OperationType::BitStore0,out,Some(nc));out
}
fn or(b: &mut B, a: BitId, c: BitId) -> BitId {
    let out=copy(b,a);bit_op(b,OperationType::BitStore1,out,Some(c));out
}
fn cx_if(b: &mut B, q: QubitId, target: QubitId, cond: BitId) {
    let mut op=Op::empty();op.kind=OperationType::CX;op.q_control1=q;
    op.q_target=target;op.c_condition=cond;b.ops.push(op);
}
fn z_if(b: &mut B, q: QubitId, cond: BitId) {
    let mut op=Op::empty();op.kind=OperationType::Z;
    op.q_target=q;op.c_condition=cond;b.ops.push(op);
}

pub(super) fn add_nbit(b: &mut B, a: &[BitId], acc: &[QubitId]) {
    let n=acc.len();assert!(a.len()<=n);if n==0{return;}
    let carries=b.alloc_qubits(n-1);
    for i in 0..n-1 {
        if i>0 {b.ccx(acc[i],carries[i-1],carries[i]);}
        if let Some(&ai)=a.get(i) {
            cx_if(b,acc[i],carries[i],ai);
            if i>0 {cx_if(b,carries[i-1],carries[i],ai);}
        }
    }
    for i in 0..n {
        if let Some(&ai)=a.get(i) {b.x_if(acc[i],ai);}
        if i>0 {b.cx(carries[i-1],acc[i]);}
    }
    for i in (0..n-1).rev() {
        let m=b.alloc_bit();b.hmr(carries[i],m);
        b.x(acc[i]);
        if i>0 {b.cz_if(acc[i],carries[i-1],m);}
        if let Some(&ai)=a.get(i) {
            let condition=and(b,ai,m);
            z_if(b,acc[i],condition);
            if i>0 {z_if(b,carries[i-1],condition);}
        }
        b.x(acc[i]);
    }
    b.free_vec(&carries);
}
pub(super) fn sub_nbit(b: &mut B, a: &[BitId], acc: &[QubitId]) {
    for &q in acc {b.x(q);}add_nbit(b,a,acc);for &q in acc {b.x(q);}
}
pub(super) fn lt_into(b: &mut B, u: &[QubitId], a: &[BitId], flag: QubitId) {
    let n=u.len();assert!(a.len()<=n);if n==0{return;}
    let carries=b.alloc_qubits(n);
    for &q in u {b.x(q);}
    for i in 0..n {
        if i>0 {b.ccx(u[i],carries[i-1],carries[i]);}
        if let Some(&ai)=a.get(i) {
            cx_if(b,u[i],carries[i],ai);
            if i>0 {cx_if(b,carries[i-1],carries[i],ai);}
        }
    }
    b.cx(carries[n-1],flag);
    for i in (0..n).rev() {
        let m=b.alloc_bit();b.hmr(carries[i],m);
        if i>0 {b.cz_if(u[i],carries[i-1],m);}
        if let Some(&ai)=a.get(i) {
            let condition=and(b,ai,m);z_if(b,u[i],condition);
            if i>0 {z_if(b,carries[i-1],condition);}
        }
    }
    for &q in u {b.x(q);}
    b.free_vec(&carries);
}
fn const_register(b: &mut B, n: usize, c: U256) -> Vec<BitId> {
    (0..n).map(|i|constant(b,i<256 && bit(c,i))).collect()
}
fn add_classical(b: &mut B, a: &[BitId], c: &[BitId], n: usize) -> Vec<BitId> {
    let zero=constant(b,false);let mut carry=zero;let mut out=Vec::with_capacity(n);
    for i in 0..n {
        let ai=a.get(i).copied().unwrap_or(zero);
        let ci=c.get(i).copied().unwrap_or(zero);
        let p=xor(b,ai,ci);out.push(xor(b,p,carry));
        let g=and(b,ai,ci);let propagated=and(b,p,carry);carry=or(b,g,propagated);
    }
    out
}
fn p_minus(b: &mut B, a: &[BitId], p: U256) -> Vec<BitId> {
    let neg:Vec<_>=a.iter().map(|&ai|not(b,ai)).collect();
    let k=const_register(b,256,p.wrapping_add(U256::from(1u64)));
    add_classical(b,&neg,&k,256)
}
pub(super) fn field_add(b: &mut B, acc: &[QubitId], a: &[BitId], p: U256) {
    assert_eq!(acc.len(),256);assert_eq!(a.len(),256);
    let (extended,overflow)=ext_reg(b,acc);
    add_nbit(b,a,&extended);
    let c=U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
    add_nbit_const_direct(b,&extended,c);
    b.x(overflow);csub_nbit_const_direct_fast(b,acc,c,overflow);b.x(overflow);
    lt_into(b,acc,a,overflow);
    unext_reg(b,overflow);
}
pub(super) fn field_sub(b: &mut B, acc: &[QubitId], a: &[BitId], p: U256) {
    assert_eq!(acc.len(),256);assert_eq!(a.len(),256);
    let threshold=p_minus(b,a,p);
    let (extended,overflow)=ext_reg(b,acc);
    sub_nbit(b,a,&extended);
    let c=U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
    csub_nbit_const_direct_fast(b,acc,c,overflow);
    b.x(overflow);lt_into(b,acc,&threshold,overflow);
    unext_reg(b,overflow);
}

pub(super) fn field_add_classical(b: &mut B, a: &[BitId], c: &[BitId], p: U256) -> Vec<BitId> {
    assert_eq!(a.len(),256);assert_eq!(c.len(),256);
    let raw=add_classical(b,a,c,257);
    let correction=U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
    let k=const_register(b,257,correction);
    let reduced=add_classical(b,&raw,&k,257);
    let flag=reduced[256];
    (0..256).map(|i| {
        let delta=xor(b,raw[i],reduced[i]);let chosen=and(b,delta,flag);xor(b,raw[i],chosen)
    }).collect()
}
