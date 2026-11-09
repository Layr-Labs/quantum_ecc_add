//! Reversible two-fold reduction of a 512-bit word modulo 2^256-c.
use super::*;

fn add_high(b: &mut B, high: &[QubitId], wide: &[QubitId], inverse: bool) {
    let terms=[(0usize,false),(4,false),(6,true),(10,false),(32,false)];
    if inverse {
        for &(shift,subtract) in terms.iter().rev() {
            clean_scale::short_term(b,high,&wide[shift..],!subtract);
        }
    } else {
        for &(shift,subtract) in &terms {
            clean_scale::short_term(b,high,&wide[shift..],subtract);
        }
    }
}

fn add_small_product(b: &mut B, q: &[QubitId], wide: &[QubitId], inverse: bool) {
    // The first fold guarantees q<=c, hence c*q<=c*c<2^65.
    // Signed intermediate terms may wrap modulo this width, but the final
    // nonnegative correction product fits without discarding information.
    let product=b.alloc_qubits(65);
    clean_scale::solinas_product(b,q,&product,false);
    if inverse {sub_zero_extended(b,&product,wide);}
    else {add_zero_extended(b,&product,wide);}
    clean_scale::solinas_product(b,q,&product,true);
    b.free_vec(&product);
}

/// Temporarily replace the low half of `word` by its canonical residue.
/// All history is retained until `body` returns, then the original word and
/// zero scratch are restored explicitly with the matching arithmetic helpers.
/// The body must preserve the word; it can use the residue as a control.
pub(super) fn with_reduced<F: FnOnce(&mut B,&[QubitId])>(
    b: &mut B, word: &[QubitId], p: U256, body: F,
) {
    assert_eq!(word.len(),512); assert_eq!(p,SECP256K1_P);
    let low=&word[..256]; let high=&word[256..];
    let c=U256::from(0x1000003d1u64);
    let q=b.alloc_qubits(33);
    let mut first=low.to_vec();first.extend_from_slice(&q);
    // T=low+c*high < 2^289. The high 33 bits are the second fold input.
    add_high(b,high,&first,false);
    let overflow=b.alloc_qubit();
    let mut sum=low.to_vec();sum.push(overflow);
    // S=T_low+c*q < 2^256+c*(c+1) < 2*p, so one reduction suffices.
    add_small_product(b,&q,&sum,false);
    let flag=b.alloc_qubit();
    add_nbit_const_direct(b,&sum,c);
    b.cx(overflow,flag);
    b.x(flag);csub_nbit_const_direct_fast(b,&sum,c,flag);b.x(flag);
    b.cx(flag,overflow);
    body(b,low);
    // Reconstruct the pre-reduction sum using the retained reduction bit.
    b.cx(flag,overflow);
    b.x(flag);cadd_nbit_const_direct_fast(b,&sum,c,flag);b.x(flag);
    b.cx(overflow,flag);
    sub_nbit_const_direct(b,&sum,c);
    b.free(flag);
    add_small_product(b,&q,&sum,true);
    b.free(overflow);
    add_high(b,high,&first,true);
    b.free_vec(&q);
}
