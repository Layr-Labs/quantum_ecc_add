use super::*;

#[test]
fn exhaustive_small_direct_arithmetic() {
    // Exhaust all pairs at small widths, packed into 64 simulation lanes.
    // This checks both arithmetic and measurement-induced phase cancellation.
    for n in 1..=6 {
        for short in [false, true] {
            if short && n == 1 { continue; }
            let aw = n - usize::from(short);
            for subtract in [false, true] {
                let mut b = B::new();
                let a = b.alloc_qubits(aw);
                let acc = b.alloc_qubits(n);
                if subtract { cuccaro_sub_fast_zero(&mut b, &a, &acc); }
                else { cuccaro_add_fast_zero(&mut b, &a, &acc); }
                check_all(&b, &a, &acc, n, subtract, None);
            }
        }
        for c in 0..1u64 << n {
            let mut b = B::new();
            let acc = b.alloc_qubits(n);
            add_nbit_const_direct(&mut b, &acc, U256::from(c));
            check_all(&b, &[], &acc, n, false, Some(c));
        }
    }
}

fn check_all(b: &B, a: &[QubitId], acc: &[QubitId], n: usize,
             subtract: bool, constant: Option<u64>) {
    let ar: Vec<_> = a.iter().copied().map(QubitOrBit::Qubit).collect();
    let rr: Vec<_> = acc.iter().copied().map(QubitOrBit::Qubit).collect();
    let cases = (1usize << a.len()) * (1usize << n);
    let mut xof = Shake256::default().finalize_xof();
    let mut sim = Simulator::new(b.next_qubit as usize, b.next_bit as usize, &mut xof);
    for first in (0..cases).step_by(64) {
        sim.clear_for_shot();
        let lanes = 64.min(cases-first);
        for lane in 0..lanes {
            let case = first+lane;
            let av = (case >> n) as u64;
            let rv = (case & ((1 << n)-1)) as u64;
            sim.set_register(&ar, U256::from(av), lane);
            sim.set_register(&rr, U256::from(rv), lane);
        }
        sim.apply_iter(b.ops.iter());
        let mask = if lanes == 64 { u64::MAX } else { (1u64 << lanes)-1 };
        assert_eq!(sim.phase & mask, 0, "phase n={n} subtract={subtract} c={constant:?}");
        for lane in 0..lanes {
            let case = first+lane;
            let av = constant.unwrap_or((case >> n) as u64);
            let rv = (case & ((1 << n)-1)) as u64;
            let expected = if subtract { rv.wrapping_sub(av) } else { rv+av } & ((1 << n)-1);
            assert_eq!(sim.get_register(&rr, lane), U256::from(expected));
            if constant.is_none() {
                assert_eq!(sim.get_register(&ar, lane), U256::from((case >> n) as u64));
            }
        }
        for q in 0..b.next_qubit {
            if !a.contains(&QubitId(q)) && !acc.contains(&QubitId(q)) {
                assert_eq!(sim.qubits[q as usize] & mask, 0, "ancilla {q}");
            }
        }
    }
}

#[test]
fn large_direct_arithmetic_and_phase() {
    let mut seed = Shake256::default();
    seed.update(b"large arithmetic coverage, no point-add test inputs");
    let mut input = seed.finalize_xof();
    for n in [8usize, 32, 255, 256] {
        let mask = if n == 256 { U256::MAX } else { (U256::from(1u64) << n) - U256::from(1u64) };
        for mode in 0..7 {
            let mut b = B::new();
            let aw = if mode == 1 || mode == 3 { n-1 } else { n };
            let a = if mode < 4 { b.alloc_qubits(aw) } else { Vec::new() };
            let acc = b.alloc_qubits(n);
            let constant = match mode {
                4 => U256::from(0x1000003d1u64),
                5 => SECP256K1_P.wrapping_add(U256::from(1u64)),
                _ => U256::MAX,
            } & mask;
            if mode < 2 { cuccaro_add_fast_zero(&mut b, &a, &acc); }
            else if mode < 4 { cuccaro_sub_fast_zero(&mut b, &a, &acc); }
            else { add_nbit_const_direct(&mut b, &acc, constant); }
            let ar: Vec<_> = a.iter().copied().map(QubitOrBit::Qubit).collect();
            let rr: Vec<_> = acc.iter().copied().map(QubitOrBit::Qubit).collect();
            let mut measurements = Shake256::default().finalize_xof();
            let mut sim = Simulator::new(b.next_qubit as usize, b.next_bit as usize, &mut measurements);
            for _ in 0..32 {
                sim.clear_for_shot();
                for &q in a.iter().chain(acc.iter()) {
                    let mut bytes = [0u8;8];
                    input.read(&mut bytes);
                    sim.qubits[q.0 as usize] = u64::from_le_bytes(bytes);
                }
                for (lane, value) in [U256::ZERO, U256::from(1u64), mask, mask-U256::from(1u64)].into_iter().enumerate() {
                    sim.set_register(&rr, value, lane);
                    let amask = if aw == 256 { U256::MAX } else { (U256::from(1u64) << aw)-U256::from(1u64) };
                    sim.set_register(&ar, value & amask, lane);
                }
                let before: Vec<_> = (0..64).map(|lane| (sim.get_register(&ar,lane),sim.get_register(&rr,lane))).collect();
                sim.apply_iter(b.ops.iter());
                assert_eq!(sim.phase,0,"large phase n={n} mode={mode}");
                for lane in 0..64 {
                    let (av,rv) = before[lane];
                    let addend = if mode<4 { av } else { constant };
                    let expected = if mode==2 || mode==3 { rv.wrapping_sub(addend) } else { rv.wrapping_add(addend) } & mask;
                    assert_eq!(sim.get_register(&rr,lane),expected,"large arithmetic n={n} mode={mode}");
                    assert_eq!(sim.get_register(&ar,lane),av);
                }
                for q in 0..b.next_qubit {
                    if !a.contains(&QubitId(q)) && !acc.contains(&QubitId(q)) { assert_eq!(sim.qubits[q as usize],0); }
                }
            }
        }
    }
}

#[test]
fn exhaustive_controlled_constants() {
    for n in 1..=5 {
        for c in 0..1u64<<n {
            for subtract in [false,true] {
                let mut b=B::new();
                let ctrl=b.alloc_qubit();
                let acc=b.alloc_qubits(n);
                if subtract { csub_nbit_const_direct_fast(&mut b,&acc,U256::from(c),ctrl); }
                else { cadd_nbit_const_direct_fast(&mut b,&acc,U256::from(c),ctrl); }
                let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
                let mut xof=Shake256::default().finalize_xof();
                let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
                sim.qubits[ctrl.0 as usize]=0xaaaaaaaaaaaaaaaa;
                for lane in 0..64 { sim.set_register(&rr,U256::from((lane>>1)&((1usize<<n)-1)),lane); }
                sim.apply_iter(b.ops.iter());
                assert_eq!(sim.phase,0,"controlled phase n={n} c={c} subtract={subtract}");
                assert_eq!(sim.qubit(ctrl),0xaaaaaaaaaaaaaaaa);
                for lane in 0..64 {
                    let rv=((lane>>1)&((1usize<<n)-1)) as u64;
                    let addend=if lane&1==1 {c} else {0};
                    let expected=if subtract {rv.wrapping_sub(addend)} else {rv+addend} & ((1u64<<n)-1);
                    assert_eq!(sim.get_register(&rr,lane),U256::from(expected));
                }
            }
        }
    }
}

#[test]
fn exhaustive_zero_extension() {
    for n in 1..=6 {
        for aw in 0..=n {
            for subtract in [false,true] {
                let mut b=B::new();
                let a=b.alloc_qubits(aw);
                let acc=b.alloc_qubits(n);
                if subtract {sub_zero_extended(&mut b,&a,&acc);}
                else {add_zero_extended(&mut b,&a,&acc);}
                check_all(&b,&a,&acc,n,subtract,None);
            }
        }
    }
}
