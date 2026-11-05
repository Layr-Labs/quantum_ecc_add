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
        apply_clean(&mut sim, &b.ops);
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
                apply_clean(&mut sim, &b.ops);
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
                apply_clean(&mut sim, &b.ops);
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

#[test]
fn exhaustive_small_karatsuba_roundtrip() {
    // The inherited full multiplier requires distinct x/y wires. Its
    // controlled add/subtract kernel is not a squaring primitive.
    for same_inputs in [false] {
        let n=4;
        let mut b=B::new();
        let x=b.alloc_qubits(n);
        let y=if same_inputs {x.clone()} else {b.alloc_qubits(n)};
        let out=b.alloc_qubits(2*n);
        let middle=b.alloc_qubits(n+2);
        karatsuba_forward(&mut b,&x,&y,&out,&middle);
        let boundary=b.ops.len();
        karatsuba_inverse(&mut b,&x,&y,&out,&middle);
        let xr:Vec<_>=x.iter().copied().map(QubitOrBit::Qubit).collect();
        let yr:Vec<_>=y.iter().copied().map(QubitOrBit::Qubit).collect();
        let rr:Vec<_>=out.iter().copied().map(QubitOrBit::Qubit).collect();
        let mut xof=Shake256::default().finalize_xof();
        let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
        for first in (0..256usize).step_by(64) {
            sim.clear_for_shot();
            for lane in 0..64 {
                let case=first+lane;
                sim.set_register(&xr,U256::from(case>>n),lane);
                if !same_inputs {sim.set_register(&yr,U256::from(case&15),lane);}
            }
            apply_clean(&mut sim, &b.ops[..boundary]);
            assert_eq!(sim.phase,0,"karatsuba phase same_inputs={same_inputs} first={first}");
            for lane in 0..64 {
                let case=first+lane;
                let xv=case>>n;
                let yv=if same_inputs {xv} else {case&15};
                assert_eq!(sim.get_register(&rr,lane),U256::from(xv*yv));
                assert_eq!(sim.get_register(&xr,lane),U256::from(xv));
                assert_eq!(sim.get_register(&yr,lane),U256::from(yv));
            }
            apply_clean(&mut sim, &b.ops[boundary..]);
            assert_eq!(sim.phase,0);
            for q in 0..b.next_qubit {
                if !x.contains(&QubitId(q)) && !y.contains(&QubitId(q)) {
                    assert_eq!(sim.qubits[q as usize],0);
                }
            }
        }
    }
}

// Reset in the simulator clears the qubit and may randomize phase. Check the
// pre-reset value explicitly so a dirty temporary cannot pass by chance.
fn apply_clean<R: XofReader>(sim: &mut Simulator<'_, R>, ops: &[Op]) {
    assert!(!ops.iter().any(|o| matches!(o.kind, OperationType::PushCondition | OperationType::PopCondition)));
    for op in ops {
        if op.kind == OperationType::R {
            assert_eq!(sim.qubit(op.q_target), 0, "dirty reset {:?}", op.q_target);
        }
        sim.apply_iter(std::iter::once(op));
    }
}

#[test]
fn exhaustive_measured_controlled_add() {
    for n in 1..=5 {
        for subtract in [false,true] {
            let mut b=B::new();
            let ctrl=b.alloc_qubit();
            let a=b.alloc_qubits(n);
            let acc=b.alloc_qubits(n);
            if subtract {controlled_sub_linear_workspace(&mut b,&a,&acc,ctrl);}
            else {controlled_add_linear_workspace(&mut b,&a,&acc,ctrl);}
            let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Qubit).collect();
            let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
            let mut xof=Shake256::default().finalize_xof();
            let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
            let cases=1usize<<(2*n+1);
            for first in (0..cases).step_by(64) {
                sim.clear_for_shot();
                let lanes=64.min(cases-first);
                for lane in 0..lanes {
                    let case=first+lane;
                    sim.qubits[ctrl.0 as usize]|=((case&1) as u64)<<lane;
                    sim.set_register(&ar,U256::from((case>>1)&((1usize<<n)-1)),lane);
                    sim.set_register(&rr,U256::from((case>>(n+1))&((1usize<<n)-1)),lane);
                }
                apply_clean(&mut sim,&b.ops);
                assert_eq!(sim.phase,0);
                for lane in 0..lanes {
                    let case=first+lane;
                    let av=((case>>1)&((1usize<<n)-1)) as u64;
                    let rv=((case>>(n+1))&((1usize<<n)-1)) as u64;
                    let addend=if case&1==1 {av} else {0};
                    let expected=if subtract {rv.wrapping_sub(addend)} else {rv+addend} & ((1u64<<n)-1);
                    assert_eq!(sim.get_register(&rr,lane),U256::from(expected),"n={n} case={case} sub={subtract}");
                    assert_eq!(sim.get_register(&ar,lane),U256::from(av));
                    assert_eq!((sim.qubit(ctrl)>>lane)&1,(case&1) as u64);
                }
            }
            assert!(b.next_qubit as usize <= 3*n);
        }
    }
}

#[test]
fn large_measured_controlled_add() {
    let mut seed=Shake256::default();
    seed.update(b"measurement controlled ladder independent large tests");
    let mut input=seed.finalize_xof();
    for occupied in [0usize,2200] {
    for n in [8usize,32,128,256] {
        for subtract in [false,true] {
            let mut b=B::new();
            let _occupied=b.alloc_qubits(occupied);
            let ctrl=b.alloc_qubit();
            let a=b.alloc_qubits(n);
            let acc=b.alloc_qubits(n);
            if subtract {controlled_sub_linear_workspace(&mut b,&a,&acc,ctrl);}
            else {controlled_add_linear_workspace(&mut b,&a,&acc,ctrl);}
            let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Qubit).collect();
            let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
            let mut xof=Shake256::default().finalize_xof();
            let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
            let mask=if n==256 {U256::MAX} else {(U256::from(1u64)<<n)-U256::from(1u64)};
            for _ in 0..32 {
                sim.clear_for_shot();
                for &q in a.iter().chain(acc.iter()).chain(std::iter::once(&ctrl)) {
                    let mut bytes=[0u8;8]; input.read(&mut bytes);
                    sim.qubits[q.0 as usize]=u64::from_le_bytes(bytes);
                }
                let cv=sim.qubit(ctrl);
                let before:Vec<_>=(0..64).map(|lane|(sim.get_register(&ar,lane),sim.get_register(&rr,lane))).collect();
                apply_clean(&mut sim,&b.ops);
                assert_eq!(sim.phase,0);
                assert_eq!(sim.qubit(ctrl),cv);
                for lane in 0..64 {
                    let (av,rv)=before[lane];
                    let addend=if (cv>>lane)&1==1 {av} else {U256::ZERO};
                    let expected=if subtract {rv.wrapping_sub(addend)} else {rv.wrapping_add(addend)} & mask;
                    assert_eq!(sim.get_register(&rr,lane),expected);
                    assert_eq!(sim.get_register(&ar,lane),av);
                }
            }
        }
    }    }
}

#[test]
fn modular_shift_measured_roundtrip() {
    let mut seed=Shake256::default(); seed.update(b"shift roundtrip independent inputs");
    let mut input=seed.finalize_xof();
    let p=SECP256K1_P;
    for k in [1usize,4,22] {
        let mut b=B::new();
        let v=b.alloc_qubits(256);
        let (spill,flag,overflow)=mod_shift_left_by_k(&mut b,&v,p,k);
        let boundary=b.ops.len();
        mod_shift_right_by_k(&mut b,&v,p,k,spill,flag,overflow);
        let vr:Vec<_>=v.iter().copied().map(QubitOrBit::Qubit).collect();
        let mut xof=Shake256::default().finalize_xof();
        let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
        for _ in 0..16 {
            sim.clear_for_shot();
            let mut before=Vec::new();
            for lane in 0..64 {
                let mut bytes=[0u8;32];input.read(&mut bytes);
                let val=match lane {
                    0=>U256::ZERO,1=>U256::from(1u64),2=>p-U256::from(1u64),
                    3=>p>>1,4=>(p>>1)+U256::from(1u64),_=>U256::from_le_bytes(bytes)%p,
                };
                before.push(val);
                sim.set_register(&vr,val,lane);
            }
            apply_clean(&mut sim,&b.ops[..boundary]);
            assert_eq!(sim.phase,0,"shift phase k={k}");
            for lane in 0..64 {
                assert_eq!(sim.get_register(&vr,lane),mulmod(before[lane],U256::from(1u64)<<k,p),"shift k={k} lane={lane}");
            }
            apply_clean(&mut sim,&b.ops[boundary..]);
            assert_eq!(sim.phase,0,"reverse shift phase k={k}");
            for lane in 0..64 {assert_eq!(sim.get_register(&vr,lane),before[lane]);}
            for q in 0..b.next_qubit {
                if !v.contains(&QubitId(q)) {assert_eq!(sim.qubits[q as usize],0);}
            }
        }
    }
}

#[test]
fn large_sparse_controlled_constants() {
    let mut seed=Shake256::default();seed.update(b"sparse controlled constant coverage");
    let mut input=seed.finalize_xof();
    for n in [8usize,32,128,256] {
        let mask=if n==256 {U256::MAX} else {(U256::from(1u64)<<n)-U256::from(1u64)};
        for c in [U256::from(1u64),U256::from(3u64),U256::from(5u64),U256::from(977u64),U256::from(0x1000003d1u64)] {
            let c=c & mask;
            for subtract in [false,true] {
                let mut b=B::new();
                let _occupied=b.alloc_qubits(2200);
                let ctrl=b.alloc_qubit();let acc=b.alloc_qubits(n);
                if subtract {csub_nbit_const_direct_fast(&mut b,&acc,c,ctrl);}
                else {cadd_nbit_const_direct_fast(&mut b,&acc,c,ctrl);}
                let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
                let mut xof=Shake256::default().finalize_xof();
                let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
                for _ in 0..16 {
                    sim.clear_for_shot();sim.qubits[ctrl.0 as usize]=0xaaaaaaaaaaaaaaaa;
                    let mut before=Vec::new();
                    for lane in 0..64 {
                        let mut bytes=[0u8;32];input.read(&mut bytes);
                        let rv=match lane/2 {0=>U256::ZERO,1=>U256::from(1u64),2=>mask,3=>c,4=>c.wrapping_sub(U256::from(1u64))&mask,_=>U256::from_le_bytes(bytes)&mask};
                        before.push(rv);sim.set_register(&rr,rv,lane);
                    }
                    apply_clean(&mut sim,&b.ops);assert_eq!(sim.phase,0);
                    assert_eq!(sim.qubit(ctrl),0xaaaaaaaaaaaaaaaa);
                    for lane in 0..64 {
                        let addend=if lane&1==1 {c} else {U256::ZERO};
                        let rv=before[lane];
                        let expected=if subtract {rv.wrapping_sub(addend)} else {rv.wrapping_add(addend)} & mask;
                        assert_eq!(sim.get_register(&rr,lane),expected,"sparse n={n} c={c} sub={subtract} lane={lane}");
                    }
                }
            }
        }
    }
}

#[test]
fn high_occupancy_field_arithmetic() {
    let p=SECP256K1_P;
    let mut seed=Shake256::default();seed.update(b"occupied field arithmetic coverage");
    let mut input=seed.finalize_xof();
    for mode in 0..3 {
        let mut b=B::new();
        let padding=b.alloc_qubits(2100);
        let a=b.alloc_qubits(256);
        let acc=b.alloc_qubits(256);
        match mode {
            0=>mod_add_qq_fast(&mut b,&acc,&a,p),
            1=>mod_sub_qq_fast(&mut b,&acc,&a,p),
            _=>mod_neg_inplace_fast(&mut b,&acc,p),
        }
        let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Qubit).collect();
        let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
        let mut xof=Shake256::default().finalize_xof();
        let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
        for _ in 0..16 {
            sim.clear_for_shot();let mut before=Vec::new();
            for lane in 0..64 {
                let mut bytes=[0u8;32];input.read(&mut bytes);
                let av=match lane/2 {0=>U256::ZERO,1=>p-U256::from(1u64),2=>U256::from(1u64),_=>U256::from_le_bytes(bytes)%p};
                input.read(&mut bytes);
                let rv=match lane/2 {0=>U256::ZERO,1=>p-U256::from(1u64),2=>p-U256::from(1u64),_=>U256::from_le_bytes(bytes)%p};
                before.push((av,rv));sim.set_register(&ar,av,lane);sim.set_register(&rr,rv,lane);
            }
            apply_clean(&mut sim,&b.ops);assert_eq!(sim.phase,0);
            for lane in 0..64 {
                let (av,rv)=before[lane];
                let expected=match mode {
                    0=>rv.add_mod(av,p),
                    1=>if rv>=av {rv-av} else {p-(av-rv)},
                    _=>p-rv,
                };
                assert_eq!(sim.get_register(&rr,lane),expected,"field mode={mode} lane={lane}");
                assert_eq!(sim.get_register(&ar,lane),av);
            }
            for q in &padding {assert_eq!(sim.qubit(*q),0);}
        }
    }
}

#[test]
fn exhaustive_classical_offset_primitives() {
    for n in 1..=6 {
        for mode in 0..3 {
            let mut b=B::new();let a=b.alloc_bits(n);let acc=b.alloc_qubits(n);let flag=b.alloc_qubit();
            match mode {0=>classical_offset::add_nbit(&mut b,&a,&acc),1=>classical_offset::sub_nbit(&mut b,&a,&acc),_=>classical_offset::lt_into(&mut b,&acc,&a,flag)}
            let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Bit).collect();
            let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
            let mut xof=Shake256::default().finalize_xof();let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
            let cases=1usize<<(2*n);
            for first in (0..cases).step_by(64) {
                sim.clear_for_shot();let lanes=64.min(cases-first);
                for lane in 0..lanes {let case=first+lane;sim.set_register(&ar,U256::from(case>>n),lane);sim.set_register(&rr,U256::from(case&((1<<n)-1)),lane);}
                apply_clean(&mut sim,&b.ops);assert_eq!(sim.phase,0);
                for lane in 0..lanes {
                    let case=first+lane;let av=(case>>n) as u64;let rv=(case&((1<<n)-1)) as u64;
                    let expected=match mode {0=>rv+av,1=>rv.wrapping_sub(av),_=>rv}&((1<<n)-1);
                    assert_eq!(sim.get_register(&rr,lane),U256::from(expected),"qb n={n} mode={mode} case={case}");
                    assert_eq!(sim.get_register(&ar,lane),U256::from(av));
                    if mode==2 {assert_eq!((sim.qubit(flag)>>lane)&1,u64::from(rv<av));}
                }
            }
        }
    }
}

#[test]
fn classical_offset_field_arithmetic() {
    let p=SECP256K1_P;let mut seed=Shake256::default();seed.update(b"classical offset field tests");let mut input=seed.finalize_xof();
    for mode in 0..3 {
        let mut b=B::new();let a=b.alloc_bits(256);let quantum=if mode<2 {b.alloc_qubits(256)}else{Vec::new()};let other=if mode==2 {b.alloc_bits(256)}else{Vec::new()};
        let result=match mode {0=>{classical_offset::field_add(&mut b,&quantum,&a,p);Vec::new()},1=>{classical_offset::field_sub(&mut b,&quantum,&a,p);Vec::new()},_=>classical_offset::field_add_classical(&mut b,&a,&other,p)};
        let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Bit).collect();
        let rr:Vec<_>=if mode<2 {quantum.iter().copied().map(QubitOrBit::Qubit).collect()}else{other.iter().copied().map(QubitOrBit::Bit).collect()};
        let out:Vec<_>=if mode<2 {rr.clone()}else{result.iter().copied().map(QubitOrBit::Bit).collect()};
        let mut xof=Shake256::default().finalize_xof();let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
        for _ in 0..16 {
            sim.clear_for_shot();let mut before=Vec::new();
            for lane in 0..64 {
                let mut bytes=[0u8;32];input.read(&mut bytes);
                let av=match lane {0=>U256::ZERO,1=>p-U256::from(1u64),2=>U256::from(1u64),_=>U256::from_le_bytes(bytes)%p};
                input.read(&mut bytes);let rv=match lane {0=>U256::ZERO,1=>p-U256::from(1u64),2=>U256::ZERO,_=>U256::from_le_bytes(bytes)%p};
                before.push((av,rv));sim.set_register(&ar,av,lane);sim.set_register(&rr,rv,lane);
            }
            apply_clean(&mut sim,&b.ops);assert_eq!(sim.phase,0);
            for lane in 0..64 {
                let (av,rv)=before[lane];let expected=if mode==1 {if rv>=av {rv-av}else{p-(av-rv)}}else{rv.add_mod(av,p)};
                assert_eq!(sim.get_register(&out,lane),expected,"qb field mode={mode} lane={lane}");assert_eq!(sim.get_register(&ar,lane),av);
                if mode==2 {assert_eq!(sim.get_register(&rr,lane),rv);}
            }
        }
    }
}

#[test]
fn classical_offset_multiples() {
    let p=SECP256K1_P;let mut b=B::new();let a=b.alloc_bits(256);
    let twice=classical_offset::field_add_classical(&mut b,&a,&a,p);
    let thrice=classical_offset::field_add_classical(&mut b,&twice,&a,p);
    assert_eq!(b.next_qubit,0);
    let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Bit).collect();
    let dr:Vec<_>=twice.iter().copied().map(QubitOrBit::Bit).collect();
    let tr:Vec<_>=thrice.iter().copied().map(QubitOrBit::Bit).collect();
    let mut seed=Shake256::default();seed.update(b"classical offset multiple coverage");let mut input=seed.finalize_xof();
    let mut xof=Shake256::default().finalize_xof();let mut sim=Simulator::new(0,b.next_bit as usize,&mut xof);
    for _ in 0..32 {
        sim.clear_for_shot();let mut before=Vec::new();
        for lane in 0..64 {
            let mut bytes=[0u8;32];input.read(&mut bytes);
            let av=match lane {0=>U256::ZERO,1=>U256::from(1u64),2=>p-U256::from(1u64),3=>p>>1,_=>U256::from_le_bytes(bytes)%p};
            before.push(av);sim.set_register(&ar,av,lane);
        }
        apply_clean(&mut sim,&b.ops);assert_eq!(sim.phase,0);assert_eq!(sim.stats.toffoli_gates,0);
        for lane in 0..64 {
            let av=before[lane];assert_eq!(sim.get_register(&ar,lane),av);
            assert_eq!(sim.get_register(&dr,lane),mulmod(av,U256::from(2u64),p));
            assert_eq!(sim.get_register(&tr,lane),mulmod(av,U256::from(3u64),p));
        }
    }
}

#[test]
fn exhaustive_constant_one_carry() {
    for n in 1..=6 {
        for short in [false,true] {
            if short && n==1 {continue;}
            for subtract in [false,true] {
                let mut b=B::new();let a=b.alloc_qubits(n-usize::from(short));let acc=b.alloc_qubits(n);
                if subtract {cuccaro_sub_fast_one(&mut b,&a,&acc);}else{cuccaro_add_fast_one(&mut b,&a,&acc);}
                let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Qubit).collect();let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
                let mut xof=Shake256::default().finalize_xof();let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
                let cases=1usize<<(a.len()+n);
                for first in (0..cases).step_by(64) {
                    sim.clear_for_shot();let lanes=64.min(cases-first);
                    for lane in 0..lanes {let case=first+lane;sim.set_register(&ar,U256::from(case>>n),lane);sim.set_register(&rr,U256::from(case&((1<<n)-1)),lane);}
                    apply_clean(&mut sim,&b.ops);assert_eq!(sim.phase,0);
                    for lane in 0..lanes {let case=first+lane;let av=(case>>n) as u64;let rv=(case&((1<<n)-1)) as u64;let expected=if subtract {rv.wrapping_sub(av+1)}else{rv+av+1}&((1<<n)-1);assert_eq!(sim.get_register(&rr,lane),U256::from(expected));assert_eq!(sim.get_register(&ar,lane),U256::from(av));}
                }
            }
        }
    }
}

#[test]
fn exhaustive_biased_controlled_kernel() {
    for n in 1..=5 {
        let mut b=B::new();let ctrl=b.alloc_qubit();let x=b.alloc_qubits(n);let acc=b.alloc_qubits(n+1);
        controlled_add_subtract_fast(&mut b,&x,&acc,ctrl);let boundary=b.ops.len();controlled_add_subtract_fast_inverse(&mut b,&x,&acc,ctrl);
        let xr:Vec<_>=x.iter().copied().map(QubitOrBit::Qubit).collect();let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
        let mut xof=Shake256::default().finalize_xof();let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
        let cases=1usize<<(2*n+2);let mask=(1u64<<(n+1))-1;
        for first in (0..cases).step_by(64) {
            sim.clear_for_shot();let lanes=64.min(cases-first);
            for lane in 0..lanes {let case=first+lane;sim.qubits[ctrl.0 as usize]|=((case&1) as u64)<<lane;sim.set_register(&xr,U256::from((case>>1)&((1<<n)-1)),lane);sim.set_register(&rr,U256::from((case>>(n+1))&mask as usize),lane);}
            apply_clean(&mut sim,&b.ops[..boundary]);assert_eq!(sim.phase,0);
            for lane in 0..lanes {let case=first+lane;let xv=((case>>1)&((1<<n)-1)) as u64;let rv=((case>>(n+1))&mask as usize) as u64;let delta=if case&1==1 {xv}else{(1u64<<n)-xv};assert_eq!(sim.get_register(&rr,lane),U256::from((rv+delta)&mask));assert_eq!(sim.get_register(&xr,lane),U256::from(xv));assert_eq!((sim.qubit(ctrl)>>lane)&1,(case&1) as u64);}
            apply_clean(&mut sim,&b.ops[boundary..]);assert_eq!(sim.phase,0);
            for lane in 0..lanes {let case=first+lane;assert_eq!(sim.get_register(&rr,lane),U256::from((case>>(n+1))&mask as usize));assert_eq!(sim.get_register(&xr,lane),U256::from((case>>1)&((1<<n)-1)));assert_eq!((sim.qubit(ctrl)>>lane)&1,(case&1) as u64);}
        }
    }
}

#[test]
fn wide_short_addend_roundtrip() {
    let mut b=B::new();let a=b.alloc_qubits(256);let acc=b.alloc_qubits(513);
    add_zero_extended(&mut b,&a,&acc);let boundary=b.ops.len();sub_zero_extended(&mut b,&a,&acc);
    let ar:Vec<_>=a.iter().copied().map(QubitOrBit::Qubit).collect();let rr:Vec<_>=acc.iter().copied().map(QubitOrBit::Qubit).collect();
    let mut seed=Shake256::default();seed.update(b"wide short addend arithmetic");let mut input=seed.finalize_xof();let mut xof=Shake256::default().finalize_xof();let mut sim=Simulator::new(b.next_qubit as usize,b.next_bit as usize,&mut xof);
    for _ in 0..16 {
        sim.clear_for_shot();
        for &q in a.iter().chain(acc.iter()) {let mut bytes=[0u8;8];input.read(&mut bytes);sim.qubits[q.0 as usize]=u64::from_le_bytes(bytes);}
        let before:Vec<_>=(0..64).map(|lane|(sim.get_register(&ar,lane),sim.get_register(&rr[..256],lane),sim.get_register(&rr[256..512],lane),(sim.qubit(acc[512])>>lane)&1)).collect();
        apply_clean(&mut sim,&b.ops[..boundary]);assert_eq!(sim.phase,0);
        for lane in 0..64 {
            let (av,lo,hi,top)=before[lane];let newlo=lo.wrapping_add(av);let carry=newlo<lo;let newhi=hi.wrapping_add(U256::from(u64::from(carry)));let highcarry=carry&&newhi<hi;
            assert_eq!(sim.get_register(&rr[..256],lane),newlo);assert_eq!(sim.get_register(&rr[256..512],lane),newhi);assert_eq!((sim.qubit(acc[512])>>lane)&1,top^u64::from(highcarry));assert_eq!(sim.get_register(&ar,lane),av);
        }
        apply_clean(&mut sim,&b.ops[boundary..]);assert_eq!(sim.phase,0);
        for lane in 0..64 {let(av,lo,hi,top)=before[lane];assert_eq!(sim.get_register(&ar,lane),av);assert_eq!(sim.get_register(&rr[..256],lane),lo);assert_eq!(sim.get_register(&rr[256..512],lane),hi);assert_eq!((sim.qubit(acc[512])>>lane)&1,top);}
    }
}
