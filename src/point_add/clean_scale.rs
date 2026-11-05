//! Clean in-place powers of two modulo the secp256k1 prime.
use super::*;

pub(super) const BLOCK_SIZE: usize = 68;

fn inverse_low_bits(c: U256, k: usize) -> U256 {
    let mut inverse=U256::from(1u64);
    for _ in 0..8 {
        inverse=inverse.wrapping_mul(U256::from(2u64).wrapping_sub(c.wrapping_mul(inverse)));
    }
    let mask=(U256::from(1u64)<<k)-U256::from(1u64);
    inverse &= mask;
    assert_eq!(inverse.wrapping_mul(c)&mask,U256::from(1u64));
    inverse
}
fn inverse_product(b: &mut B, value: &[QubitId], q: &[QubitId], inverse: U256, subtract: bool) {
    let k=q.len();let mask=(U256::from(1u64)<<k)-U256::from(1u64);
    for i in 0..k {
        let constant=(inverse<<i)&mask;
        if subtract {csub_nbit_const_direct_fast(b,q,constant,value[i]);}
        else {cadd_nbit_const_direct_fast(b,q,constant,value[i]);}
    }
}
fn short_term(b: &mut B, q: &[QubitId], target: &[QubitId], subtract: bool) {
    assert!(target.len()>=q.len());
    let padding=b.alloc_qubits(target.len()-q.len());
    let mut input=q.to_vec();input.extend_from_slice(&padding);
    if subtract {sub_nbit_qq_fast(b,&input,target);}else{add_nbit_qq_fast(b,&input,target);}
    b.free_vec(&padding);
}
fn solinas_product(b: &mut B, q: &[QubitId], product: &[QubitId], inverse: bool) {
    // c = 1 + 2^4 - 2^6 + 2^10 + 2^32. Products fit k+33 bits.
    let terms=[(4usize,false),(6,true),(10,false),(32,false)];
    if !inverse {
        for i in 0..q.len() {b.cx(q[i],product[i]);}
        for &(shift,negative) in &terms {short_term(b,q,&product[shift..],negative);}
    } else {
        for &(shift,negative) in terms.iter().rev() {short_term(b,q,&product[shift..],!negative);}
        for i in 0..q.len() {b.cx(q[i],product[i]);}
    }
}
fn below_product_into(b: &mut B, value: &[QubitId], product: &[QubitId], target: QubitId) {
    let w=product.len();assert!(w<value.len());
    let high_zero=b.alloc_qubit();let low_less=b.alloc_qubit();
    with_eq_zero_fast(b,&value[w..],high_zero,|b| {
        with_lt(b,&value[..w],product,low_less,|b| {
            b.ccx(high_zero,low_less,target);
        });
    });
    b.free(low_less);b.free(high_zero);
}

pub(super) fn left(b: &mut B, value: &[QubitId], p: U256, k: usize) {
    assert_eq!(p,SECP256K1_P);assert_eq!(value.len(),256);assert!(k>0&&k<=128);
    let c=U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
    let inverse=inverse_low_bits(c,k);
    let (q,flag,overflow)=mod_shift_left_by_k(b,value,p,k);
    // flag is !reduce; the old top-bit spill plus reduce is the quotient q.
    b.x(flag);cadd_nbit_const_direct_fast(b,&q,U256::from(1u64),flag);
    let product=b.alloc_qubits(k+33);
    solinas_product(b,&q,&product,false);
    below_product_into(b,value,&product,flag);
    solinas_product(b,&q,&product,true);b.free_vec(&product);
    b.free(flag);b.free(overflow);
    // value mod 2^k = c*q mod 2^k, so subtract its known inverse product.
    inverse_product(b,value,&q,inverse,true);b.free_vec(&q);
}

pub(super) fn right(b: &mut B, value: &[QubitId], p: U256, k: usize) {
    assert_eq!(p,SECP256K1_P);assert_eq!(value.len(),256);assert!(k>0&&k<=128);
    let c=U256::MAX.wrapping_sub(p).wrapping_add(U256::from(1u64));
    let inverse=inverse_low_bits(c,k);
    let q=b.alloc_qubits(k);inverse_product(b,value,&q,inverse,false);
    let flag=b.alloc_qubit();let product=b.alloc_qubits(k+33);
    solinas_product(b,&q,&product,false);
    below_product_into(b,value,&product,flag);
    solinas_product(b,&q,&product,true);b.free_vec(&product);
    // Recover the original captured top bits and the inverse reduction flag.
    csub_nbit_const_direct_fast(b,&q,U256::from(1u64),flag);b.x(flag);
    let overflow=b.alloc_qubit();
    mod_shift_right_by_k(b,value,p,k,q,flag,overflow);
}

pub(super) fn many(b: &mut B, value: &[QubitId], p: U256, steps: usize, inverse: bool) {
    let mut remaining=steps;
    while remaining>0 {
        let k=remaining.min(BLOCK_SIZE);
        if inverse {right(b,value,p,k);}else{left(b,value,p,k);}
        remaining-=k;
    }
}
