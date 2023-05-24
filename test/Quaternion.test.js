import Quaternion from "../src/Quaternion";
import AMath from "../src/AMath";
import Vector3 from "../src/Vector3";

test('should work with different params', () => {
    let quat = new Quaternion();

    expect(quat.w).toEqual(1);
    expect((quat.add(quat)).w).toEqual(2);
    expect(new Quaternion(2).w).toEqual(2);
    expect(new Quaternion(2).x).toEqual(0);

    expect(new Quaternion(2, 3).w).toEqual(2);
    expect(new Quaternion(2, 3).x).toEqual(3);
    expect(new Quaternion(2, 3).y).toEqual(0);
    expect(new Quaternion(2, 3).z).toEqual(0);
})

test("should work with arrays", () => {
    expect((new Quaternion(1, [2, 3, 4])).w).toEqual(1);
    expect((new Quaternion(1, [2, 3, 4])).x).toEqual(2);
    expect((new Quaternion(1, [2, 3, 4])).y).toEqual(3);
    expect((new Quaternion(1, [2, 3, 4])).z).toEqual(4);

    expect((new Quaternion([1, 2, 3, 4])).w).toEqual(1);
    expect((new Quaternion([1, 2, 3, 4])).x).toEqual(2);
    expect((new Quaternion([1, 2, 3, 4])).y).toEqual(3);
    expect((new Quaternion([1, 2, 3, 4])).z).toEqual(4);

    expect((new Quaternion([2, 3, 4])).w).toEqual(2);
    expect((new Quaternion([2, 3, 4])).x).toEqual(3);
    expect((new Quaternion([2, 3, 4])).y).toEqual(4);
    expect((new Quaternion([2, 3, 4])).z).toEqual(0);
})

test("should parse strings", () => {
    expect(new Quaternion(1).toString()).toEqual('1');
    expect(new Quaternion(1 + 1).toString()).toEqual('2');
    expect(new Quaternion(1 - 1).toString()).toEqual('0');
})

test('should return correct vectors', () => {
    expect(Quaternion.zero.toString()).toEqual('0');
    expect(Quaternion.one.toString()).toEqual('1');
    expect(Quaternion.i.toString()).toEqual('i');
    expect(Quaternion.j.toString()).toEqual('j');
    expect(Quaternion.k.toString()).toEqual('k');
    expect(Quaternion.epsilon.toString()).toEqual('1e-16');
})

test('should convert from string to quaternion', () => {
    expect((Quaternion.fromString('1, 4, 5, 8')).w).toEqual(1);
    expect((Quaternion.fromString('1, 4, 5, 8')).x).toEqual(4);
    expect((Quaternion.fromString('1, 4, 5, 8')).y).toEqual(5);
    expect((Quaternion.fromString('1, 4, 5, 8')).z).toEqual(8);

    expect((Quaternion.fromString('1+2i-4j+1k')).w).toEqual(1);
    expect((Quaternion.fromString('1+2i-4j+1k')).x).toEqual(2);
    expect((Quaternion.fromString('1+2i-4j+1k')).y).toEqual(-4);
    expect((Quaternion.fromString('1+2i-4j+1k')).z).toEqual(1);
    expect((Quaternion.fromString('1+2i-4j+k')).y).toEqual(-4);
    expect((Quaternion.fromString('1+2i-4j+k')).z).toEqual(1);
})

test('should add two quats', () => {
    expect(((new Quaternion(1, 2, 3, 4)).add(new Quaternion(9, 8, 7, -6))).toString()).toEqual('10 + 10i + 10j - 2k');
    expect(((new Quaternion(1, -2, 3, -4)).add(new Quaternion(9, 8, 7, 6))).toString()).toEqual('10 + 6i + 10j + 2k');
    expect(((new Quaternion(-9, -2, 3, -4)).add(new Quaternion(9, 8, 7, 6))).toString()).toEqual('6i + 10j + 2k');
    expect(new Quaternion().add(new Quaternion()).toString()).toEqual('2');
    expect((new Quaternion(2)).add(new Quaternion(-1)).toString()).toEqual('1');

    expect((new Quaternion(1)).add(new Quaternion(0, 1))).toEqual(new Quaternion(1, 1, 0, 0))
    
    expect((new Quaternion(1, 2, 3, 4)).add(new Quaternion(5, 6, 7, 8)).toString()).toEqual('6 + 8i + 10j + 12k');
    expect((new Quaternion(-1, 2, 3, 4)).add(new Quaternion(5, 6, 7, 8)).toString()).toEqual('4 + 8i + 10j + 12k');
    expect((new Quaternion(1, -2, 3, 4)).add(new Quaternion(5, 6, 7, 8)).toString()).toEqual('6 + 4i + 10j + 12k');
    expect((new Quaternion(1, 2, -3, 4)).add(new Quaternion(5, 6, 7, 8)).toString()).toEqual('6 + 8i + 4j + 12k');
    expect((new Quaternion(1, 2, 3, -4)).add(new Quaternion(5, 6, 7, 8)).toString()).toEqual('6 + 8i + 10j + 4k');

    expect(Quaternion.zero.add(Quaternion.zero)).toEqual(new Quaternion());
    expect((new Quaternion(1)).add(new Quaternion(-1))).toEqual(Quaternion.zero);
    expect((new Quaternion(0, 0, 1)).add(Quaternion.zero)).toEqual(new Quaternion(0, 0, 1, 0));
    expect((new Quaternion(0, 0, 0, 1)).add(Quaternion.zero)).toEqual(new Quaternion(0, 0, 0, 1));
})

test("should subtract two quats", () => {
    expect(((new Quaternion(1, 2, 3, 4)).sub(new Quaternion(9, 8, 4, -6))).toString()).toEqual('-8 - 6i - j + 10k');
    expect(((new Quaternion(1, -2, 3, -4)).sub(new Quaternion(9, 8, 4, 6))).toString()).toEqual('-8 - 10i - j - 10k');
    expect(((new Quaternion(-9, -2, 3, -4)).sub(new Quaternion(9, 8, 4, 6))).toString()).toEqual('-18 - 10i - j - 10k');
    expect((Quaternion.zero.sub(Quaternion.zero)).toString()).toEqual('0');
    expect(((new Quaternion(2)).sub(new Quaternion(-1))).toString()).toEqual('3');
    expect(((new Quaternion(0, 0, 0, 1)).sub(new Quaternion(-1))).toString()).toEqual('1 + k');

    expect(Quaternion.zero.sub(Quaternion.zero)).toEqual(Quaternion.zero);
    expect(Quaternion.zero.sub(new Quaternion(1, 2, 3, 4))).toEqual(new Quaternion(-1, -2, -3, -4));
    expect((new Quaternion(10, 9, 8, 7)).sub(new Quaternion(1, 2, 3, 4))).toEqual(new Quaternion(9, 7, 5, 3));
})

test("should calculate the norm of a quat", () => {
    expect((new Quaternion()).norm()).toEqual(1);
    expect((new Quaternion(1, 1, 1, 1)).norm()).toEqual(2);
    expect((new Quaternion([3, 2, 5, 4])).normalize().norm()).toEqual(1);
    expect((new Quaternion(1, 2, 3, 4)).norm()).toEqual(Math.sqrt(1 + 4 + 9 + 16));
    expect((new Quaternion(1, 2, 3, 4)).normSq()).toEqual(1 + 4 + 9 + 16);

    expect((new Quaternion(5)).norm()).toEqual(5);
    expect((new Quaternion(-5)).norm()).toEqual(5);
    expect((new Quaternion(1, 1, 1, 1)).norm()).toEqual(2);

    expect((new Quaternion(0)).norm()).toEqual(0);
    expect((new Quaternion(3, 4)).norm()).toEqual(5);
    expect((new Quaternion(-3, 4)).norm()).toEqual(5);
    expect((new Quaternion(0, -3, 4)).norm()).toEqual(5);
    expect((new Quaternion(0, 0, -3, 4)).norm()).toEqual(5);
    expect((new Quaternion(1, 1, 1, 1)).norm()).toEqual(2);

    expect((new Quaternion(1, 2, 6, 20)).norm()).toEqual(21);
    expect((new Quaternion(20, 1, 2, 6)).norm()).toEqual(21);

    // expect((new Quaternion('i')).fromString().norm()).toEqual(1)
    expect((new Quaternion([3, 2, 5, 4])).norm()).toBeCloseTo(7.348);
})

test('should calculate the inverse of a quat', () => {
    expect((new Quaternion(1, 2, 3, 4)).inverse().toString()).toBe('0.03333333333333333 - 0.06666666666666667i - 0.1j - 0.13333333333333333k');

    let p = new Quaternion([3, 2, 5, 4]);
    let p_ = p.conjugate();
    let l = p.norm();
    let r = 1 / (l * l);

    expect(p_.norm()).toBeCloseTo(l);
    expect(p_.scale(r)).toEqual(p.inverse());
})

test('should calculate the conjugate of a quat', () => {
    expect((new Quaternion(1, 2, 3, 4)).conjugate().toString()).toEqual('1 - 2i - 3j - 4k');

    expect((new Quaternion(1, 2, 3, 4)).conjugate().w).toEqual(1);
    expect((new Quaternion(1, 2, 3, 4)).conjugate().x).toEqual(-2);
    expect((new Quaternion(1, 2, 3, 4)).conjugate().y).toEqual(-3);
    expect((new Quaternion(1, 2, 3, 4)).conjugate().z).toEqual(-4);

    expect((new Quaternion(0, 0, 0, 0)).conjugate().w).toEqual(0);
    expect((new Quaternion(0, 0, 0, 0)).conjugate().x).toEqual(-0);
    expect((new Quaternion(0, 0, 0, 0)).conjugate().y).toEqual(-0);
    expect((new Quaternion(0, 0, 0, 0)).conjugate().z).toEqual(-0);

    expect((new Quaternion(1)).conjugate()).toEqual(new Quaternion([1, -0, -0, -0]));
    expect((new Quaternion('i')).conjugate()).toEqual(new Quaternion([0, -1, -0, -0]));
    expect((new Quaternion('j')).conjugate()).toEqual(new Quaternion([0, -1, -0, -0]));
    expect((new Quaternion('k')).conjugate()).toEqual(new Quaternion([0, -1, -0, -0]));
})

test('should pass conjugate properties', () => {
    let p1 = new Quaternion(8, 1, 2, 3);
    let p2 = new Quaternion(6, 9, 8, 7);

    expect(p2.conjugate().mul(p1.conjugate())).toEqual(p1.mul(p2).conjugate());

    let p = new Quaternion(Math.random(), Math.random(), Math.random(), Math.random()).normalize();

    expect(p.mul(p.conjugate())).toEqual(p.conjugate().mul(p));

    var a = new Quaternion(1);

    expect(a.conjugate()).toEqual(a);

    var a = new Quaternion(0);

    expect(a.conjugate()).toEqual(a);

    let q1 = new Quaternion(0, 1, 2, 3);
    let q2 = new Quaternion(0, 9, 8, 7);

    expect(q1.mul(q2).conjugate()).toEqual(q2.mul(q1));
})

test('should pass hamilton rules', () => {
    let i2 = (new Quaternion('i')).mul(new Quaternion('i'));
    let j2 = (new Quaternion('j')).mul(new Quaternion('j'));
    let k2 = (new Quaternion('k')).mul(new Quaternion('k'));
    let ijk = (new Quaternion('i')).mul(new Quaternion('j')).mul(new Quaternion('k'));

    expect(j2).toEqual(i2);
    expect(k2).toEqual(j2);
    expect(ijk).toEqual(k2);
    expect(new Quaternion(-1, 0, 0, 0)).toEqual(ijk);

    let qI = new Quaternion(0, 1);
    let qJ = new Quaternion(0, 0, 1);
    let qK = new Quaternion(0, 0, 0, 1);

    expect(qI.mul(qI)).toEqual(new Quaternion(-1));
    expect(qJ.mul(qJ)).toEqual(new Quaternion(-1));
    expect(qK.mul(qK)).toEqual(new Quaternion(-1));

    expect(qI.mul(qJ)).toEqual(new Quaternion("k"));
    expect(qJ.mul(qI)).toEqual(new Quaternion("-k"));
    expect(qJ.mul(qK)).toEqual(new Quaternion("i"));
    expect(qK.mul(qJ)).toEqual(new Quaternion("-i"));
    expect(qK.mul(qI)).toEqual(new Quaternion("j"));
    expect(qI.mul(qK)).toEqual(new Quaternion("-j"));
})

test('should add a number to a Quaternion', () => {
    expect((new Quaternion(1, 2, 3, 4)).add(5)).toEqual(new Quaternion(6, 2, 3, 4));
    expect((new Quaternion(1, 2, 3, 4)).add(-5)).toEqual(new Quaternion(-4, 2, 3, 4));
    expect((new Quaternion(1, 2, 3, 4)).add(0)).toEqual(new Quaternion(1, 2, 3, 4));
    expect((new Quaternion(0, 0, 0, 0)).add(5)).toEqual(new Quaternion(5, 0, 0, 0));
})

test('should return the real and imaginary part', () => {
    let q = new Quaternion(7, 2, 3, 4);

    expect((new Quaternion(q.imag())).toString()).toEqual('2i + 3j + 4k');
    expect(q.real()).toEqual(7);
})

test('should result in the same for the inverse of normalized quats', () => {
    let q = new Quaternion(9, 8, 7, 6).normalize();

    expect(q.inverse()).toEqual(q.conjugate());
})

test('should normalize quaternion', () => {
    let q = (new Quaternion(Math.random() * 1000, Math.random() * 1000, Math.random() * 1000, Math.random() * 1000)).normalize();

    expect(new Quaternion(q.norm())).toEqual(new Quaternion(1, 0, 0, 0));
})

test('should invert quaternion', () => {
    let q = new Quaternion(Math.random() * 1000, Math.random() * 1000, Math.random() * 1000, Math.random() * 1000);

    expect(q.mul(q.inverse())).toEqual(new Quaternion(1, 0, 0, 0));
    expect(q.inverse().mul(q)).toEqual(new Quaternion(1, 0, 0, 0));
})

test('should calculate the dot product', () => {
    expect((new Quaternion(9, 8, 7, 6)).dot(new Quaternion(1, 2, 3, 4)).toString()).toEqual('70');
    expect((new Quaternion(9, 8, 7, 6)).normSq()).toEqual((new Quaternion(9, 8, 7, 6)).dot(new Quaternion(9, 8, 7, 6)));
})

test('should pass trivial cases', () => {
    var q0 = new Quaternion(0);
    var q1 = new Quaternion(Math.random(), Math.random(), Math.random(), Math.random());
    var q2 = new Quaternion(Math.random(), Math.random(), Math.random(), Math.random());
    var l = Math.random() * 2.0 - 1.0;
    var lp = Math.random();

    expect(q1.add(q2)).toEqual(q2.add(q1));
    expect(q0.sub(q1)).toEqual(q1.neg());
    expect(q1.conjugate().conjugate()).toEqual(q1);
    expect(q1.normalize().norm()).toBeCloseTo(1);
    expect(q1.inverse()).toEqual(q1.conjugate().scale(1 / Math.pow(q1.norm(), 2)));
    expect(q1.div(q2)).toEqual(q1.mul(q2.inverse()));
    expect(q1.mul(q2).norm()).toBeCloseTo(q1.norm() * q2.norm());
    expect((new Quaternion(l)).exp()).toEqual(new Quaternion(Math.exp(l)));
    expect((new Quaternion(lp)).log()).toEqual(new Quaternion(Math.log(lp)));
    expect(q1.exp().log()).toEqual(q1);
    expect(q1.log().exp()).toEqual(q1);
    expect(q1.pow(2.0)).toEqual(q1.mul(q1));
    expect(q1.mul(q1.inverse())).toEqual(Quaternion.one);
    expect(q1.inverse().mul(q1)).toEqual(Quaternion.one);
    expect(q1.add(q1.conjugate())).toEqual(new Quaternion(2 * q1.w));
})

test('should pass other trivial cases', () => {
    var x = new Quaternion(1, 2, -0.5, -1);
    var y = new Quaternion(-3, 4, 0, 2);
    var z = new Quaternion(-2, 1, 2, -4);

    expect(y.normSq()).toEqual(29);
    expect(z.normSq()).toEqual(25);
    expect(z.normalize()).toEqual(new Quaternion(-0.4, 0.2, 0.4, -0.8));
    expect(x.exp().log()).toEqual(x);
    expect(x.mul(y)).toEqual(new Quaternion(-9.0, -3.0, -6.5, 7.0));
    expect(y.dot(y)).toEqual(29);
})

test('should calculate the product', () => {
    expect(new Quaternion(5).mul(new Quaternion(6)).toString()).toEqual('30');
    expect(new Quaternion(1, 2, 3, 4).mul(new Quaternion(6)).toString()).toEqual('6 + 12i + 18j + 24k'); // scale
    expect(new Quaternion(6).mul(new Quaternion(1, 2, 3, 4)).toString()).toEqual('6 + 12i + 18j + 24k'); // scale
    expect(new Quaternion(5, 6).mul(new Quaternion(6, 7)).toString()).toEqual('-12 + 71i');
    expect(new Quaternion(1, 1, 1, 1).mul(new Quaternion(2, 2, 2, 2)).toString()).toEqual('-4 + 4i + 4j + 4k');

    expect(new Quaternion(1, 2, 3, 4).mul(new Quaternion(5, 6, 7, 8))).toEqual(new Quaternion(-60, 12, 30, 24));
    expect(new Quaternion(3, 2, 5, 4).mul(new Quaternion(4, 5, 3, 1))).toEqual(new Quaternion(-17, 16, 47, 0));
    expect(new Quaternion().mul(new Quaternion(1, 2, 3, 4))).toEqual(new Quaternion(1, 2, 3, 4));
    expect(new Quaternion().mul(new Quaternion())).toEqual(new Quaternion());

    expect(new Quaternion(1, 0, 0, 0).mul(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(1, 0, 0, 0));
    expect(new Quaternion(0, 1, 0, 0).mul(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(0, 1, 0, 0));
    expect(new Quaternion(0, 0, 1, 0).mul(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(0, 0, 1, 0));
    expect(new Quaternion(0, 0, 0, 1).mul(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(0, 0, 0, 1));
    expect(new Quaternion(1, 0, 0, 0).mul(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(0, 1, 0, 0));
    expect(new Quaternion(0, 1, 0, 0).mul(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(-1, 0, 0, 0));
    expect(new Quaternion(0, 0, 1, 0).mul(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(0, 0, 0, -1));
    expect(new Quaternion(0, 0, 0, 1).mul(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(0, 0, 1, 0));
    expect(new Quaternion(1, 0, 0, 0).mul(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(0, 0, 1, 0));
    expect(new Quaternion(0, 1, 0, 0).mul(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(0, 0, 0, 1));
    expect(new Quaternion(0, 0, 1, 0).mul(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(-1, 0, 0, 0));
    expect(new Quaternion(0, 0, 0, 1).mul(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(0, -1, 0, 0));
})

test('should scale a quaternion', () => {
    expect((new Quaternion([3, 2, 5, 4])).scale(3)).toEqual(new Quaternion([9, 6, 15, 12]));
})

test('should calculate the quotient', () => {
    expect((new Quaternion(6)).div(new Quaternion(2)).toString()).toEqual('3');
    expect((new Quaternion(1)).div(new Quaternion(2)).toString()).toEqual('0.5');
    expect((new Quaternion(1)).div(new Quaternion(2)).toString()).toEqual('0.5');
    expect((new Quaternion(4, 2)).div(new Quaternion(1, 1)).toString()).toEqual('3 - i');
    expect((new Quaternion(3, -2)).div(Quaternion.i).toString()).toEqual('-2 - 3i');
    expect((new Quaternion(25)).div(new Quaternion(3, -4)).toString()).toEqual('3 + 4i');
})

test('should result in norm=1 with fromAxisAngle', () => {
    var axis = [1, 1, 1];
    var angle = Math.PI;

    expect(Quaternion.fromAxisAngle(axis, angle).norm()).toEqual(1);
})

test('should have no effect to rotate on axis parallel to vector direction', () => {
    var v = [1, 1, 1];

    var angle = Math.random();
    var axis = [1, 1, 1];

    var r = Quaternion.fromAxisAngle(axis, angle).rotateVector(v);

    expect(r).toEqual([1,1,1]);
})

test('should generate a unit quaternion from euler angle', () => {
    var n = Quaternion.fromEuler(Math.PI, Math.PI, Math.PI).norm();

    expect(n).toEqual(1);
})

test("should rotate a vector in direct and indirect manner", () => {
    var v = [1, 9, 3];

    var q = (new Quaternion("1+2i+3j+4k")).normalize();

    var a = q.mul(v).mul(q.conjugate()).toVector();
    var b = q.rotateVector(v);

    expect(a.slice(1)).toEqual(b);
});

test("should rotate a vector correctly", () => {
    var theta = 2 * Math.PI / 3.0;
    var axis = [1.0, 1.0, 1.0];
    var vector = [3.0, 4.0, 5.0];

    var p = Quaternion.fromAxisAngle(axis, theta).rotateVector(vector);

    expect(p).toEqual([5.0, 3.0, 4.0]);
});

test("should rotate a vector correctly", () => {
    var v = [1.0, 1.0, 1.0];
    var q = Quaternion.fromAxisAngle([0.0, 1.0, 0.0], Math.PI);
    var p = q.rotateVector(v);

    expect(p[0]).toBeCloseTo(-1);
    expect(p[1]).toBeCloseTo(1);
    expect(p[2]).toBeCloseTo(-1);
});

test("should rotate a vector correctly based on Euler angles I", () => {
    var v = [1.0, 2.0, 3.0];

    var q = Quaternion.fromEuler(0.0, Math.PI, 0.0, 'ZXY');
    var p = q.rotateVector(v);
    expect(p[0]).toBeCloseTo(1);
    expect(p[1]).toBeCloseTo(-2);
    expect(p[2]).toBeCloseTo(-3);

    var q = Quaternion.fromEuler(Math.PI, 0.0, 0.0, 'ZXY');
    var p = q.rotateVector(v);
    expect(p[0]).toBeCloseTo(-1);
    expect(p[1]).toBeCloseTo(-2);
    expect(p[2]).toBeCloseTo(3);

    var q = Quaternion.fromEuler(Math.PI, 0, Math.PI/2, 'ZXY');
    var p = q.rotateVector(v);
    expect(p[0]).toBeCloseTo(-3);
    expect(p[1]).toBeCloseTo(-2);
    expect(p[2]).toBeCloseTo(-1);
});

test("should rotate a vector correctly based on Euler angles II", () => {
    var v = [1.0, 2.0, 3.0];

    var q = Quaternion.fromEuler(0.0, Math.PI, 0.0, 'XYZ');
    var p = q.rotateVector(v);

    expect(p[0]).toBeCloseTo(-1);
    expect(p[1]).toBeCloseTo(2);
    expect(p[2]).toBeCloseTo(-3);

    var q = Quaternion.fromEuler(Math.PI, 0.0, 0.0, 'XYZ');
    var p = q.rotateVector(v);

    expect(p[0]).toBeCloseTo(1);
    expect(p[1]).toBeCloseTo(-2);
    expect(p[2]).toBeCloseTo(-3);

    var q = Quaternion.fromEuler(Math.PI, 0, Math.PI/2, 'XYZ');
    var p = q.rotateVector(v);

    expect(p[0]).toBeCloseTo(-2);
    expect(p[1]).toBeCloseTo(-1);
    expect(p[2]).toBeCloseTo(-3);
});

test("should rotate a vector correctly based on Euler angles II", () => {
    var v = [1.0, 2.0, 3.0];

    var q = Quaternion.fromEuler(0.0, Math.PI, 0.0, 'XZY');
    var p = q.rotateVector(v);

    expect(p[0]).toBeCloseTo(-1);
    expect(p[1]).toBeCloseTo(-2);
    expect(p[2]).toBeCloseTo(3);

    var q = Quaternion.fromEuler(Math.PI, 0.0, 0.0, 'XZY');
    var p = q.rotateVector(v);

    expect(p[0]).toBeCloseTo(1);
    expect(p[1]).toBeCloseTo(-2);
    expect(p[2]).toBeCloseTo(-3);

    var q = Quaternion.fromEuler(Math.PI, 0, Math.PI/2, 'XZY');
    var p = q.rotateVector(v);

    expect(p[0]).toBeCloseTo(3);
    expect(p[1]).toBeCloseTo(-2);
    expect(p[2]).toBeCloseTo(1);
});

test("should exp and log a quaternion", () => {
    var q = new Quaternion(Math.random() * 10, Math.random() * 10, Math.random() * 10, Math.random() * 10);

    expect(q.log().exp()).toEqual(q);
});

test("should exp and log real numbers", () => {
    var n = Math.random() * 10;
    var q = new Quaternion(n);

    expect(q.exp().toVector()).toEqual([Math.exp(n), 0, 0, 0]);
    expect(q.log().toVector()).toEqual([Math.log(n), 0, 0, 0]);
});

test("should work with scalar powers", () => {
    var q = new Quaternion(Math.random() * 10, Math.random() * 10, Math.random() * 10, Math.random() * 10);

    expect(q.mul(q).mul(q)).toEqual(q.pow(3));
});

test('should square Quaternions', () => {
    expect((new Quaternion("i")).pow(2)).toEqual(new Quaternion(-1));
    expect((new Quaternion("j")).pow(2)).toEqual(new Quaternion(-1));
    expect((new Quaternion("k")).pow(2)).toEqual(new Quaternion(-1));
    expect((new Quaternion(1)).pow(2)).toEqual(new Quaternion({w: 1}));

    expect((new Quaternion(3, -2, -3, 4)).pow(2)).toEqual(new Quaternion(-20, -12, -18, 24));
    expect((new Quaternion(1, 2, 3, 4)).pow(2)).toEqual(new Quaternion(-28, 4, 6, 8));
    expect((new Quaternion(-1, 2, 3, 4)).pow(2)).toEqual(new Quaternion(-28, -4, -6, -8));
    expect((new Quaternion(1, -2, 3, 4)).pow(2)).toEqual(new Quaternion(-28, -4, 6, 8));
    expect((new Quaternion(1, 2, -3, 4)).pow(2)).toEqual(new Quaternion(-28, 4, -6, 8));
    expect((new Quaternion(1, 2, 3, -4)).pow(2)).toEqual(new Quaternion(-28, 4, 6, -8));

    expect((new Quaternion(5, 4, 3, 2)).pow(2)).toEqual(new Quaternion(-4, 40, 30, 20));
    expect((new Quaternion(-5, 4, 3, 2)).pow(2)).toEqual(new Quaternion(-4, -40, -30, -20));
    expect((new Quaternion(5, -4, 3, 2)).pow(2)).toEqual(new Quaternion(-4, -40, 30, 20));
    expect((new Quaternion(5, 4, -3, 2)).pow(2)).toEqual(new Quaternion(-4, 40, -30, 20));
    expect((new Quaternion(5, 4, 3, -2)).pow(2)).toEqual(new Quaternion(-4, 40, 30, -20));
});

test('should raise Quaternion to a Quaternion power', () => {
    expect((new Quaternion(1, 4, 0, 0)).pow(new Quaternion(-2, 3, 0, 0))).toEqual(new Quaternion(-0.000030177061938851806, 0.0011015451057806702, 0, 0));
    expect((new Quaternion(1, 4, -3, 2)).pow(new Quaternion(4, 2, -3, 2))).toEqual(new Quaternion(4.023822744421112, -0.08808649248602358, 0.10799947333843203, -0.045858528052467734));
    expect((new Quaternion(-1, -1, 0, 4)).pow(new Quaternion(-4, 5, 1, 1.5))).toEqual(new Quaternion(0.00009562614911354535, 0.0010196374737841477, 0.0015348157881126755, -0.0007464390363321687));

    let q1 = (new Quaternion(0, 2, 0, 0)).pow(new Quaternion(1, 0, 0, 0));
    let q2 = new Quaternion(0, 2, 0, 0);

    expect(q1.w).toEqual(q2.w);
    expect(q1.x).toEqual(q2.x);
    expect(q1.y).toEqual(q2.y);
    expect(q1.z).toEqual(q2.z);
    expect((new Quaternion(0, 2, 0, 0)).pow(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(0.15990905692806803, 0.13282699942462048, 0, 0));
    expect((new Quaternion(0, 2, 0, 0)).pow(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(-0.145615694487965, 0, 0.399409670603132, 0.905134235650981));
    expect((new Quaternion(0, 2, 0, 0)).pow(new Quaternion(0, 0, 0, 1))).toEqual(new Quaternion(-0.145615694487965, 0, -0.905134235650981, 0.399409670603132));

    expect((new Quaternion(0, 0, 2, 0)).pow(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(0, 0, 2, 0));
    expect((new Quaternion(0, 0, 2, 0)).pow(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(-0.145615694487965, 0.399409670603132, 0, -0.905134235650981));
    expect((new Quaternion(0, 0, 2, 0)).pow(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(0.159909056928068, 0, 0.13282699942462, 0));
    expect((new Quaternion(0, 0, 2, 0)).pow(new Quaternion(0, 0, 0, 1))).toEqual(new Quaternion(-0.145615694487965, 0.905134235650981, 0, 0.399409670603132));

    expect((new Quaternion(0, 0, 0, 2)).pow(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(0, 0, 0, 2));
    expect((new Quaternion(0, 0, 0, 2)).pow(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(-0.145615694487965, 0.399409670603132, 0.905134235650981, 0));
    expect((new Quaternion(0, 0, 0, 2)).pow(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(-0.145615694487965, -0.905134235650981, 0.399409670603132, 0));
    expect((new Quaternion(0, 0, 0, 2)).pow(new Quaternion(0, 0, 0, 1))).toEqual(new Quaternion(0.159909056928068, 0, 0, 0.13282699942462));
});

test('should raise reals to quaternion powers', () => {
    expect((new Quaternion(1)).pow(new Quaternion(3, 4, 5, 9))).toEqual(new Quaternion(1));
    expect((new Quaternion(-2)).pow(new Quaternion(4, 2, 1, 1.5))).toEqual(new Quaternion(-0.024695944127665907, 0.015530441791896946, -0.004473740387907085, 0.004654139181719533));
    expect((new Quaternion(2, 0, 0, 0)).pow(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(2, 0, 0, 0));
    expect((new Quaternion(2, 0, 0, 0)).pow(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(0.7692389013639721, 0.6389612763136348, 0, 0));
    expect((new Quaternion(2, 0, 0, 0)).pow(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(0.769238901363972, 0, 0.638961276313635, 0));
    expect((new Quaternion(2, 0, 0, 0)).pow(new Quaternion(0, 0, 0, 1))).toEqual(new Quaternion(0.769238901363972, 0, 0, 0.638961276313635));
});

test('should return the square root of a Quaternion', () => {
    expect((new Quaternion(1, 2, 3, 4)).pow(new Quaternion(1 / 2))).toEqual(new Quaternion(1.7996146219471072, 0.5556745248702426, 0.8335117873053639, 1.1113490497404852));
    expect((new Quaternion(-1, -2, -3, -4)).pow(new Quaternion(1 / 2)).pow(new Quaternion(2))).toEqual(new Quaternion(-1, -2, -3, -4));
    expect((new Quaternion()).pow(new Quaternion(1 / 2))).toEqual(new Quaternion());
    expect((new Quaternion(1, 0, 0, 0)).pow(new Quaternion(1 / 2))).toEqual(new Quaternion(1, 0, 0, 0));
    expect((new Quaternion(0, 1, 0, 0)).pow(new Quaternion(1 / 2))).toEqual(new Quaternion(0.7071067811865476, 0.7071067811865475, 0, 0));
    expect((new Quaternion("1-2i")).pow(new Quaternion(1 / 2))).toEqual(new Quaternion("1.272019649514069 - 0.7861513777574233i"))
    expect((new Quaternion(-1)).pow(new Quaternion(1 / 2))).toEqual(new Quaternion("i"));
});

test('should return the log base e of a quaternion', () => {
    expect((new Quaternion(1, 2, 3, 4)).log()).toEqual(new Quaternion(1.7005986908310777, 0.515190292664085, 0.7727854389961275, 1.03038058532817));
    expect((new Quaternion(-1, 2, 3, 4)).log()).toEqual(new Quaternion(1.7005986908310777, 0.6515679277817118, 0.9773518916725678, 1.3031358555634236));
    expect((new Quaternion(1, -2, 3, 4)).log()).toEqual(new Quaternion(1.7005986908310777, -0.515190292664085, 0.7727854389961275, 1.03038058532817));
    expect((new Quaternion(1, 2, -3, 4)).log()).toEqual(new Quaternion(1.7005986908310777, 0.515190292664085, -0.7727854389961275, 1.03038058532817));
    expect((new Quaternion(1, 2, 3, -4)).log()).toEqual(new Quaternion(1.7005986908310777, 0.515190292664085, 0.7727854389961275, -1.03038058532817));

    expect((new Quaternion(2, 3, 4, 5)).log()).toEqual(new Quaternion(1.9944920232821373, 0.549487105217117, 0.7326494736228226, 0.9158118420285283));
    expect((new Quaternion(5, 2, 3, 4)).log()).toEqual(new Quaternion(1.9944920232821373, 0.30545737557546476, 0.45818606336319717, 0.6109147511509295));
    expect((new Quaternion(4, 5, 2, 3)).log()).toEqual(new Quaternion(1.9944920232821373, 0.8072177296195943, 0.3228870918478377, 0.48433063777175656));
    expect((new Quaternion(3, 4, 5, 2)).log()).toEqual(new Quaternion(1.9944920232821373, 0.685883734654061, 0.8573546683175763, 0.3429418673270305));
});

test('should return the exp of a quaternion', () => {
    expect((new Quaternion(0, 0, 0, 0)).exp()).toEqual(new Quaternion(1, 0, 0, 0));
    expect((new Quaternion(1, 0, 0, 0)).exp()).toEqual(new Quaternion(2.718281828459045, 0, 0, 0));
    expect((new Quaternion(0, 1, 0, 0)).exp()).toEqual(new Quaternion(0.5403023058681398, 0.8414709848078965, 0, 0));
    expect((new Quaternion(0, 0, 1, 0)).exp()).toEqual(new Quaternion(0.5403023058681398, 0, 0.8414709848078965, 0));
    expect((new Quaternion(0, 0, 0, 1)).exp()).toEqual(new Quaternion(0.5403023058681398, 0, 0, 0.8414709848078965));

    expect((new Quaternion(-1, 0, 0, 0)).exp()).toEqual(new Quaternion(0.3678794411714424, 0, 0, 0));
    expect((new Quaternion(0, -1, 0, 0)).exp()).toEqual(new Quaternion(0.5403023058681398, -0.8414709848078965, 0, 0));
    expect((new Quaternion(0, 0, -1, 0)).exp()).toEqual(new Quaternion(0.5403023058681398, 0, -0.8414709848078965, 0));
    expect((new Quaternion(0, 0, 0, -1)).exp()).toEqual(new Quaternion(0.5403023058681398, 0, 0, -0.8414709848078965));

    expect((new Quaternion(1, 2, 3, 4)).exp()).toEqual(new Quaternion(1.6939227236832994, -0.7895596245415588, -1.184339436812338, -1.5791192490831176));
    expect((new Quaternion(4, 1, 2, 3)).exp()).toEqual(new Quaternion(-45.05980201339819, -8.240025266756877, -16.480050533513754, -24.720075800270628));
    expect((new Quaternion(3, 4, 1, 2)).exp()).toEqual(new Quaternion(-2.6000526954284027, -17.384580348249628, -4.346145087062407, -8.692290174124814));
    expect((new Quaternion(2, 3, 4, 1)).exp()).toEqual(new Quaternion(2.786189997492657, -4.026439818820405, -5.3685864250938735, -1.3421466062734684));
});

test('should divide quaternions by each other', () => {
    expect((new Quaternion({z: 1})).div(new Quaternion({y: 1}))).toEqual(new Quaternion({x: 1}));
    expect((new Quaternion({x: 1})).div(new Quaternion({z: 1}))).toEqual(new Quaternion({y: 1}));
    expect((new Quaternion(3, -2, -3, 4)).div(new Quaternion(3, -2, -3, 4))).toEqual(new Quaternion(1, 0, 0, 0));
    expect((new Quaternion(1, 2, 3, 4)).div(new Quaternion(-1, 1, 2, 3))).toEqual(new Quaternion(19 / 15, -4 / 15, -1 / 5, -8 / 15));

    expect((new Quaternion(1, 0, 0, 0)).div(new Quaternion(1, 0, 0, 0))).toEqual(new Quaternion(1, 0, 0, 0));
    expect((new Quaternion(1, 0, 0, 0)).div(new Quaternion(0, 1, 0, 0))).toEqual(new Quaternion(0, -1, 0, 0));
    expect((new Quaternion(1, 0, 0, 0)).div(new Quaternion(0, 0, 1, 0))).toEqual(new Quaternion(0, 0, -1, 0));
    expect((new Quaternion(1, 0, 0, 0)).div(new Quaternion(0, 0, 0, 1))).toEqual(new Quaternion(0, 0, 0, -1));
});

test('should raise Quaternion to real powers', () => {
    expect((new Quaternion(1, 2, 3, 4)).pow(2)).toEqual(new Quaternion(-28, 4, 6, 8));
    expect((new Quaternion(1, 2, 3, 4)).pow(0)).toEqual(new Quaternion({w: 1}));
    expect((new Quaternion(1, 2, 3, 4)).pow(2.5)).toEqual(new Quaternion(-66.50377063575604, -8.360428208578368, -12.54064231286755, -16.720856417156735));
    expect((new Quaternion(1, 2, 3, 4)).pow(-2.5)).toEqual(new Quaternion(-0.0134909686430938, 0.0016959981926818065, 0.0025439972890227095, 0.003391996385363613));
});

test('should rotate the optimized function', () => {
    var v = [Math.random() * 100, Math.random() * 50, Math.random() * 20];
    var q = Quaternion.random();
    let quat1 = q.mul(new Quaternion(0, v)).mul(q.conjugate()).toVector().slice(1);
    let quat2 = q.rotateVector(v);

    expect(quat1[0]).toBeCloseTo(quat2[0]);
    expect(quat1[1]).toBeCloseTo(quat2[1]);
    expect(quat1[2]).toBeCloseTo(quat2[2]);
});

test('should rotate one vector onto the other', () => {
    var u = [Math.random() * 100, Math.random() * 100, Math.random() * 100];
    var v = [Math.random() * 100, Math.random() * 100, Math.random() * 100];

    var q = Quaternion.fromBetweenVectors(u, v);
    var vPrime = q.rotateVector(u);

    // Is the length of rotated equal to the original?
    expect((new Quaternion(u)).norm()).toBeCloseTo((new Quaternion(vPrime)).norm());

    // Do they look in the same direction?
    expect((new Quaternion(v)).normalize()).toEqual((new Quaternion(vPrime)).normalize());
});

test('should rotate additive inverse to the same point', () => {
    var q1 = new Quaternion(Math.random(), Math.random(), Math.random(), Math.random()).normalize();

    var v = [Math.random(),Math.random(),Math.random()];

    expect(q1.neg().rotateVector(v)).toEqual(q1.rotateVector(v));
    expect(q1.conjugate().neg().rotateVector(v)).toEqual(q1.conjugate().rotateVector(v));
});

test('should slerp around', () => {
    var q1 = new Quaternion(Math.random(),Math.random(),Math.random(),Math.random()).normalize();
    var q2 = new Quaternion(Math.random(),Math.random(),Math.random(),Math.random()).normalize();

    expect(q1.slerp(q2)(0)).toEqual(q1);
    expect(q1.slerp(q2)(1)).toEqual(q2);
});

test('Should create a quaternion from a matrix with m22 positive and m00 < -m11', () => {
    const q = Quaternion.fromMatrix([
        -0.14040120936120104, -0.03817338250185204,   0.9893585261563572,
        -0.2992237727316569,  -0.9508943214366381,    -0.0791525318090902, 
        0.9437969242597403,   -0.3071527019707341,    0.12208432917426992 ])
    expect(q).toEqual(new Quaternion(-0.08773368562933877, 0.6496939246485931, -0.12982927130494584, 0.7438716051799714))
})

test('Should create a quaternion from a matrix with m22 positive and m00 >= -m11', () => {
    const q = Quaternion.fromMatrix([
        0.7441469599075261, -0.6679403460655629,  0.010049684482756005,
        0.4879070091702578, 0.5331745589946937,   -0.6911379312722953,
        0.4562806729009248, 0.5192115019321414,   0.7226530037289329 ])
    expect(q).toEqual(new Quaternion(0.8660217264351906, 0.34939926916920544, -0.12881633762671046, 0.3336658076679094))
})

test('Should create a matrix from a quaternion with m22 positive and m00 >= -m11', () => {
    let matrix1 = [[-0.1332, -0.0358,  0.9854],
        [-0.3022, -0.9402, -0.0754],
        [0.9386, -0.3094, 0.1212]];

    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[0][0]).toBeCloseTo(matrix1[0][0]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[0][1]).toBeCloseTo(matrix1[0][1]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[0][2]).toBeCloseTo(matrix1[0][2]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[1][0]).toBeCloseTo(matrix1[1][0]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[1][1]).toBeCloseTo(matrix1[1][1]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[1][2]).toBeCloseTo(matrix1[1][2]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[2][0]).toBeCloseTo(matrix1[2][0]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[2][1]).toBeCloseTo(matrix1[2][1]);
    expect((new Quaternion(-0.09, 0.65, -0.13, 0.74)).toMatrix(true)[2][2]).toBeCloseTo(matrix1[2][2]);
})

test('Should create a matrix from a quaternion with m22 negative and m00 > m11', () => {
    let matrix1 = [[-0.1404, -0.0382, 0.9894], 
        [-0.2992, -0.9509, -0.0792], 
        [0.9438, -0.3072, 0.1221]];
    
    let input = (new Quaternion(-0.08773368562933877, 0.6496939246485931, -0.12982927130494584, 0.7438716051799714)).toMatrix(true);
    expect(input[0][0]).toBeCloseTo(matrix1[0][0]);
    expect(input[0][1]).toBeCloseTo(matrix1[0][1]);
    expect(input[0][2]).toBeCloseTo(matrix1[0][2]);
    expect(input[1][0]).toBeCloseTo(matrix1[1][0]);
    expect(input[1][1]).toBeCloseTo(matrix1[1][1]);
    expect(input[1][2]).toBeCloseTo(matrix1[1][2]);
    expect(input[2][0]).toBeCloseTo(matrix1[2][0]);
    expect(input[2][1]).toBeCloseTo(matrix1[2][1]);
    expect(input[2][2]).toBeCloseTo(matrix1[2][2]);
})

it('should convert from and to matrix', function() {
    var initialQuat = Quaternion.random().normalize();

    let m1 = initialQuat.toMatrix(true);
    let m2 = initialQuat.toMatrix(false);

    let q1 = Quaternion.fromMatrix(m1)
    let q2 = Quaternion.fromMatrix(m2)

    expect((q1.toMatrix(true))[0][0]).toBeCloseTo[m1[0][0]];
    expect((q1.toMatrix(true))[0][1]).toBeCloseTo[m1[0][1]];
    expect((q1.toMatrix(true))[0][2]).toBeCloseTo[m1[0][2]];
    expect((q1.toMatrix(true))[1][0]).toBeCloseTo[m1[1][0]];
    expect((q1.toMatrix(true))[1][1]).toBeCloseTo[m1[1][1]];
    expect((q1.toMatrix(true))[1][2]).toBeCloseTo[m1[1][2]];
    expect((q1.toMatrix(true))[2][0]).toBeCloseTo[m1[2][0]];
    expect((q1.toMatrix(true))[2][1]).toBeCloseTo[m1[2][1]];
    expect((q1.toMatrix(true))[2][2]).toBeCloseTo[m1[2][2]];

    expect((q2.toMatrix(false))[0]).toBeCloseTo[m2[0]];
    expect((q2.toMatrix(false))[1]).toBeCloseTo[m2[1]];
    expect((q2.toMatrix(false))[2]).toBeCloseTo[m2[2]];
    expect((q2.toMatrix(false))[3]).toBeCloseTo[m2[3]];
    expect((q2.toMatrix(false))[4]).toBeCloseTo[m2[4]];
    expect((q2.toMatrix(false))[5]).toBeCloseTo[m2[5]];
    expect((q2.toMatrix(false))[6]).toBeCloseTo[m2[6]];
    expect((q2.toMatrix(false))[7]).toBeCloseTo[m2[7]];
    expect((q2.toMatrix(false))[8]).toBeCloseTo[m2[8]];
});

// fromEuler tests

test('Should fuzz 0# ZXY', () => {
    let expectExp = Quaternion.fromEuler(0.7560671181764893, -0.13546265539440805, 3.0647447738212223, 'ZXY').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-2.73473153490426) 
    expect(expectExp[1]).toBeCloseTo(0.552994141668246) 
    expect(expectExp[2]).toBeCloseTo(-5.49685736683069) 
});

test('Should fuzz 1# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-0.9380873336181326, -1.5456730703056945, -3.0597760267130223, 'ZXY').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-5.24518381954551);
    expect(expectExp[1]).toBeCloseTo(-0.867660165532979);
    expect(expectExp[2]).toBeCloseTo(-3.12013021143754);
});

test('Should fuzz 2# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-1.1348507711853153, -0.4115641368061649, 0.9705915513821521, 'ZXY').rotateVector([2, 3, 5]); 

    expect(expectExp[0]).toBeCloseTo(5.13724048032114);
    expect(expectExp[1]).toBeCloseTo(-3.40488714700346);
    expect(expectExp[2]).toBeCloseTo(-0.124514109724319);
});

test('Should fuzz 3# ZXY', () => {
    let expectExp = Quaternion.fromEuler(1.4089894888325434, 0.8042682209887531, -0.7023692389520657, 'ZXY').rotateVector([2, 3, 5]);
    
    expect(expectExp[0]).toBeCloseTo(1.30362097726652); 
    expect(expectExp[1]).toBeCloseTo(-1.93885362106741); 
    expect(expectExp[2]).toBeCloseTo(5.70450865401259);
});

test('Should fuzz 4# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-0.0344911398695098, -0.6950067407079987, -2.4174145744101403, 'ZXY').rotateVector([2, 3, 5]);
    
    expect(expectExp[0]).toBeCloseTo(-4.78181622048923);
    expect(expectExp[1]).toBeCloseTo(0.919731535018256);
    expect(expectExp[2]).toBeCloseTo(-3.77999041492953); 
});

test('Should fuzz 5# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-0.3993324281430666, -0.8686914211519521, -2.1443171860341117, 'ZXY').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-4.42266701817088); 
    expect(expectExp[1]).toBeCloseTo(3.11332137523277);
    expect(expectExp[2]).toBeCloseTo(-2.95757442187043); 
});

test('Should fuzz 6# ZXY', () => {
    let expectExp = Quaternion.fromEuler(0.7842113974343352, 1.729805218326823, -2.180977197459205, 'ZXY').rotateVector([2, 3, 5]);
    
    expect(expectExp[0]).toBeCloseTo(-4.23175395401849);
    expect(expectExp[1]).toBeCloseTo(-3.18279350186644); 
    expect(expectExp[2]).toBeCloseTo(3.15627692022193); 
});

test('Should fuzz 7# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-0.028137277945027517, 2.559355041887808, 0.7286830611984181, 'ZXY').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(4.71203393180756);
    expect(expectExp[1]).toBeCloseTo(-3.95874898470177); 
    expect(expectExp[2]).toBeCloseTo(-0.353613774642442); 
});

test('Should fuzz 8# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-3.0680638615270226, 0.24400957228642461, -2.5937427905734074, 'ZXY').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(4.57103196205473);
    expect(expectExp[1]).toBeCloseTo(-3.36393470288412);
    expect(expectExp[2]).toBeCloseTo(-2.40616086673481); 
});

test('Should fuzz 9# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-2.6316946431369574, -1.9757245652714845, -0.4350628450289924, 'ZXY').rotateVector([2, 3, 5]);
    
    expect(expectExp[0]).toBeCloseTo(2.09172605948023);
    expect(expectExp[1]).toBeCloseTo(-3.13877002139904);
    expect(expectExp[2]).toBeCloseTo(-4.87573633873469);
});

test('Should fuzz 10# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-0.058132427889624694, 1.6446132037760925, -0.6053880104987908, 'ZXY').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-1.51583698019482);
    expect(expectExp[1]).toBeCloseTo(-5.37753704154801); 
    expect(expectExp[2]).toBeCloseTo(2.60467533797456); 
});

test('Should fuzz 11# ZXY', () => {
    let expectExp = Quaternion.fromEuler(-1.3717044853059401, 2.346223855318967, -1.625263027398635, 'ZXY').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-4.27495589642530);
    expect(expectExp[1]).toBeCloseTo(4.34173147193671);
    expect(expectExp[2]).toBeCloseTo(0.934943800029096);
});


test('Should fuzz 12# XYZ', () => {
    let expectExp = Quaternion.fromEuler(2.710678164003398, 2.205556479654592, 2.054694398057457, 'XYZ').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(6.15253010525126);
    expect(expectExp[1]).toBeCloseTo(-0.308095113846526);
    expect(expectExp[2]).toBeCloseTo(0.226827478055163);
});

test('Should fuzz 13# XYZ', () => {
    let expectExp = Quaternion.fromEuler(-1.3496927969636376, -1.1026224650513043, -0.0594490015216973, 'XYZ').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-3.48061570987999);
    expect(expectExp[1]).toBeCloseTo(4.72550979783212);
    expect(expectExp[2]).toBeCloseTo(-1.88543666844826);
});

test('Should fuzz 14# XYZ', () => {
    let expectExp = Quaternion.fromEuler(-0.46335342031213766, 3.0524058758259427, -0.4280831633095632, 'XYZ').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-2.60738496003823);
    expect(expectExp[1]).toBeCloseTo(-0.649084452072119);
    expect(expectExp[2]).toBeCloseTo(-5.54799360528153);
});

test('Should fuzz 15# XYZ', () => {
    let expectExp = Quaternion.fromEuler(-0.548287375026244, -2.6546411123852907, -0.7230893888395862, 'XYZ').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-5.41926941538940);
    expect(expectExp[1]).toBeCloseTo(-0.663132616199000);
    expect(expectExp[2]).toBeCloseTo(-2.86212755424319);
});


test('Should fuzz 16# YXZ', () => {
    let expectExp = Quaternion.fromEuler(-2.02305178889919, -1.1994191611651226, -1.8724281406460914, 'YXZ').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-4.97182595580561); 
    expect(expectExp[1]).toBeCloseTo(3.64268300172091);
    expect(expectExp[2]).toBeCloseTo(0.108661005660800); 
});

test('Should fuzz 17# YXZ', () => {
    let expectExp = Quaternion.fromEuler(-2.6569753209313918, 0.10265235953947549, 2.4981777112146712, 'YXZ').rotateVector([2, 3, 5]);
    
    expect(expectExp[0]).toBeCloseTo(0.748623129852677);
    expect(expectExp[1]).toBeCloseTo(-1.70633532092011);
    expect(expectExp[2]).toBeCloseTo(-5.87605166604498);
});

test('Should fuzz 18# YXZ', () => {
    let expectExp = Quaternion.fromEuler(0.27232821692530695, -2.6879786200457056, 2.8317312488047346, 'YXZ').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(-3.65960463109698);
    expect(expectExp[1]).toBeCloseTo(4.21109035513079);
    expect(expectExp[2]).toBeCloseTo(-2.62183370276950);
});

test('Should fuzz 19# YXZ', () => {
    let expectExp = Quaternion.fromEuler(-2.573089664363896, 2.087031268358495, 1.5930548937937559, 'YXZ').rotateVector([2, 3, 5]);

    expect(expectExp[0]).toBeCloseTo(2.98880328621701);
    expect(expectExp[1]).toBeCloseTo(-5.30243457650877);
    expect(expectExp[2]).toBeCloseTo(-0.975316604053605);
});