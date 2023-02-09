import Vector3 from "../src/Vector3";

test('check getting coordinates from vector', () => {
    let vec0 = Vector3.zero;
    let vec1 = new Vector3(2, 3, 4);

    expect(vec0.x).toBe(0);
    expect(vec0.y).toBe(0);
    expect(vec0.z).toBe(0);
    expect(vec1.x).toBe(2);
    expect(vec1.y).toBe(3);
    expect(vec1.z).toBe(4);
})

test('default vectors', () => {
    expect(Vector3.zero).toEqual(new Vector3(0, 0, 0));
    expect(Vector3.right).toEqual(new Vector3(1, 0, 0));
    expect(Vector3.up).toEqual(new Vector3(0, 1, 0));
    expect(Vector3.forward).toEqual(new Vector3(0, 0, 1));
})

test('check magnitude and magnitude squared calculation of the vector', () => {
    let vec1 = new Vector3(1, 1, 1);
    let vec2 = new Vector3(0, 0, 0);
    let vec3 = new Vector3(-2, -3, -3);
    let vec4 = new Vector3(Math.sqrt(9), Math.sqrt(2), Math.sqrt(5));

    expect(vec1.magnitude).toEqual(Math.sqrt(3));
    expect(vec1.magnitudeSquared).toEqual(3);
    expect(vec2.magnitude).toEqual(0);
    expect(vec2.magnitudeSquared).toEqual(0);
    expect(vec3.magnitude).toEqual(Math.sqrt(22));
    expect(vec3.magnitudeSquared).toEqual(22);
    expect(vec4.magnitude).toEqual(4);
    expect(vec4.magnitudeSquared).toEqual(16);
})

test('check normalize property', () => {
    let vec1 = new Vector3(1, 1, 1);
    let vec2 = new Vector3(0, 0, 0);
    let vec3 = new Vector3(-2, -3, -3);
    let vec4 = new Vector3(7, 8, 6);

    expect(vec2.normalize()).toEqual(new Vector3(null, null, null));
    expect(vec1.normalize()).toEqual(new Vector3(1 / Math.sqrt(3), 1 / Math.sqrt(3), 1 / Math.sqrt(3)))
    expect(vec3.normalize()).toEqual(
        new Vector3(vec3.x / vec3.magnitude, vec3.y / vec3.magnitude, vec3.z / vec3.magnitude)
    )
    expect(vec4.normalize()).toEqual(
        new Vector3(vec4.x / vec4.magnitude, vec4.y / vec4.magnitude, vec4.z / vec4.magnitude)
    )
})

test('check neg property', () => {
    expect(new Vector3(9, 2, -1).neg()).toEqual(new Vector3(-9, -2, 1));
    expect(new Vector3(0, 0, 0).neg()).toEqual(new Vector3(0, 0, 0));
})

test('check add function', () => {
    let vec0 = Vector3.zero;
    let vec1 = new Vector3(-2, 1, 3);
    let vec2 = new Vector3(-4, 5, 8);
    let vecSum = Vector3.add(new Vector3(1, -1, 1), new Vector3(-5, 4, -3));
    expect(Vector3.add(vec0, vec0)).toEqual(vec0);
    expect(Vector3.add(vec1, vec2)).toEqual(new Vector3(-6, 6, 11));
    expect(Vector3.add(vecSum, vec1)).toEqual(new Vector3(-8, 8, 6));
})

test('check sub function', () => {
    let vec0 = Vector3.zero;
    let vec1 = new Vector3(-5, -9, -8);
    let vec2 = new Vector3(-4, 5, 8);
    let vecSum = Vector3.add(vec1, vec2); //-4;3;-2
    expect(Vector3.sub(vec1, vec1)).toEqual(Vector3.zero);
    expect(Vector3.sub(vec0, vec2)).toEqual(new Vector3(4, -5, -8));
    expect(Vector3.sub(vecSum, vec1)).toEqual(vec2);
})

test('check mul function', () => {
    let vec0 = Vector3.zero;
    let vec1 = new Vector3(4, 3, 2);
    let vec2 = new Vector3(-2, -3, -4);
    let vec3 = new Vector3(-4, 5, -6);
    expect(Vector3.mul(vec0, Math.random())).toEqual(Vector3.zero);
    expect(Vector3.mul(vec1, 0)).toEqual(Vector3.zero);
    expect(Vector3.mul(vec1, -1)).toEqual(new Vector3(-4, -3, -2));
    expect(Vector3.mul(vec2, -1)).toEqual(new Vector3(2, 3, 4));
    expect(Vector3.mul(vec3, 3)).toEqual(new Vector3(-12, 15, -18));
    expect(Vector3.mul(vec2, -10)).toEqual(new Vector3(20, 30, 40));
})

test('check div function', () => {
    let vec0 = Vector3.zero;
    let vec1 = new Vector3(4, 3, 2);
    let vec2 = new Vector3(-2, -3, -4);
    let vec3 = new Vector3(-4, 5, -6);
    expect(Vector3.div(vec0, 0)).toEqual(new Vector3(null, null, null));
    expect(Vector3.div(vec1, 0)).toEqual(new Vector3(null, null, null));
    expect(Vector3.div(vec0, 2)).toEqual(vec0);
    expect(Vector3.div(vec2, -2)).toEqual(new Vector3(1, 1.5, 2));
    expect(Vector3.div(vec3, -20)).toEqual(new Vector3(0.2, -0.25, 0.3));
    expect(Vector3.div(vec3, 20)).toEqual(new Vector3(-0.2, 0.25, -0.3));
})

test('check dot product function', () => {
    let vec0 = new Vector3(0, 2, 0);

    expect(Vector3.dot(vec0, new Vector3(1, 0, 0))).toEqual(0);
    expect(Vector3.dot(new Vector3(-6, -7, 8), new Vector3(2, -3, -4))).toEqual(-23);
    expect(Vector3.dot(vec0, new Vector3(9, 1, -7))).toEqual(2);
})

test('check cross product function', () => {
    let vec0 = new Vector3(0, 2, 0);
    let vec1 = new Vector3(1, 0, 0);
    let vec3 = new Vector3(2, -3, 3);
    let vec4 = new Vector3(4, -6, 6);
    let vec5 = new Vector3(1, -8, -7);

    expect(Vector3.cross(vec0, vec1)).toEqual(new Vector3(0, 0, -2));
    expect(Vector3.cross(vec3, vec4)).toEqual(Vector3.zero);
    expect(Vector3.cross(vec3, vec5)).toEqual(new Vector3(45, 17, -13));
    expect(Vector3.cross(vec4, vec5)).toEqual(new Vector3(90, 34, -26));
})

test('check angleBetween function', () => {
    let vec0 = new Vector3(0, 2, 0);
    let vec1 = new Vector3(1, 0, 0);
    let vec2 = new Vector3(2, -3, 3);
    let vec3 = new Vector3(4, -6, 6);
    let vec4 = new Vector3(1, -8, -7);

    expect(Vector3.angleBetween(vec0, vec1)).toEqual(Math.PI / 2);
    expect(Vector3.angleBetween(vec2, vec3)).toEqual(0);
    expect(Vector3.angleBetween(vec2, vec4)).toBeCloseTo(
        Math.acos(5 / (Math.sqrt(22) * Math.sqrt(114))),
        1);
    expect(Vector3.angleBetween(vec3, vec4)).toBeCloseTo(
        Math.acos(10 / (Math.sqrt(76) * Math.sqrt(114))),
        1);
})

test('check angleBetweenSigned function', () => {
    let vec0 = new Vector3(0, 2, 0);
    let vec1 = new Vector3(1, 0, 0);
    let vec2 = new Vector3(2, -3, 3);
    let vec3 = new Vector3(4, -6, 6);
    let vec4 = new Vector3(1, -8, -7);

    expect(Vector3.angleBetweenSigned(vec0, vec1, new Vector3(0, 1, 0))).toEqual(0);
    expect(Vector3.angleBetweenSigned(vec0, vec1, new Vector3(0, 0, 1))).toBeCloseTo(Vector3.angleBetween(vec0, vec1) * (-1));
    expect(Vector3.angleBetweenSigned(vec2, vec3, new Vector3(1, 0, 0))).toEqual(0);
    expect(Vector3.angleBetweenSigned(vec2, vec4, new Vector3(1, 0, 0))).toBeCloseTo(
        Math.acos(5 / (Math.sqrt(22) * Math.sqrt(114))),
        1);
    expect(Vector3.angleBetweenSigned(vec2, vec4, new Vector3(0, 1, 0))).toBeCloseTo(
        Math.acos(5 / (Math.sqrt(22) * Math.sqrt(114))),
        1);
    expect(Vector3.angleBetweenSigned(vec2, vec4, new Vector3(0, 0, 1))).toBeCloseTo(
        Math.acos(5 / (Math.sqrt(22) * Math.sqrt(114))) * (-1),
        1);
    expect(Vector3.angleBetweenSigned(vec2, vec4, new Vector3(0, 0, 0))).toEqual(0);
})

test('check distance function', () => {
    let vec0 = new Vector3(1, 2, 3);
    let vec1 = new Vector3(-4, -6, 0);
    let vec2 = new Vector3(9, -8, 7);

    expect(typeof Vector3.distance(vec0, vec1)).toEqual('number');
    expect(Vector3.distance(vec0, vec0)).toEqual(0);
    expect(Vector3.distance(vec0, vec1)).toEqual(Math.sqrt(98));
    expect(Vector3.distance(vec0, vec2)).toEqual(Math.sqrt(180));
    expect(Vector3.distance(vec1, vec2)).toEqual(Math.sqrt(222));
})

test('check vector rotation function', () => {
    expect(Vector3.rotate(new Vector3(1, 2, 3), { w: 0.5, x: 2, y: 4, z: 6 })).toEqual(new Vector3(42.25, 84.5, 122.75));
    expect(Vector3.rotate(new Vector3(0, 0, 0), { w: 2, x: 2, y: 1, z: -2 })).toEqual(new Vector3(0, 0, 0));
    expect(Vector3.rotate(new Vector3(-2, -2, -3), { w: 0, x: 1, y: 2, z: 0 })).toEqual(new Vector3(-22, -34, -15));
})

test('check projectOnPlane function', () => {
    expect(Vector3.projectOnPlane(new Vector3(2, 1, 1), new Vector3(0, 1, 0))).toEqual(new Vector3(2, 1, 1));
    expect(Vector3.projectOnPlane(new Vector3(2, 1, 1), new Vector3(1, 0, 0))).toEqual(new Vector3(4, 2, 2));
    expect(Vector3.projectOnPlane(new Vector3(2, 1, 1), new Vector3(0, 0, 1))).toEqual(new Vector3(2, 1, 1));
    expect(Vector3.projectOnPlane(new Vector3(2, 1, 1), new Vector3(0, 0, 0))).toEqual(new Vector3(2, 1, 1));
    expect(Vector3.projectOnPlane(new Vector3(-1, -3, -5), new Vector3(-2, 4, 0))).toEqual(new Vector3(0.5, 1.5, 2.5));
    expect(Vector3.projectOnPlane(new Vector3(2, -5, -8), new Vector3(1, 1, 1))).toEqual(new Vector3(-2, 5, 8));
})

test('check linearly interpolation between two vectors (lerp function)', () => {
    expect(Vector3.lerp(new Vector3(-2, 4, 2), new Vector3(1, -5, 2), 0.6)).toEqual(new Vector3(-0.2, -1.4, 2));
    expect(Vector3.lerp(new Vector3(-2, 4, 2), new Vector3(1, -5, 2), -0.6)).toEqual(new Vector3(-2, 4, 2));
    expect(Vector3.lerp(new Vector3(-2, 4, 2), new Vector3(1, -5, 2), 2)).toEqual(new Vector3(1, -5, 2));
    expect(Vector3.lerp(new Vector3(0, 1, 0), new Vector3(1, 0, 0), 0.5)).toEqual(new Vector3(0.5, 0.5, 0.5));
})

test('check clamMagnitude property', () => {
    let vec0 = new Vector3(0, 0, 0);
    let vec1 = new Vector3(0, 4, 0);
    let vec2 = new Vector3(Math.sqrt(10), Math.sqrt(6), 3);

    expect(vec0.clampMagnitude(-9)).toEqual(new Vector3(null, null, null));
    expect(vec0.clampMagnitude(1)).toEqual(new Vector3(0, 0, 0));
    expect(vec1.clampMagnitude(1)).toEqual(new Vector3(0, 1, 0));
    expect(vec2.clampMagnitude(2.5)).toEqual(
        new Vector3(Math.sqrt(10) / 2, Math.sqrt(6) / 2, 1.5)
    )
})

test('check fromArray method', () => {
    expect(Vector3.fromArray([0, 2, 3]) instanceof Vector3).toBeTruthy();
    expect(Vector3.fromArray([0, 2, 3])).toEqual(new Vector3(0, 2, 3));
})

test('check toString property', () => {
    expect(typeof new Vector3(8, 9, 10).toString()).toEqual('string');
    expect(new Vector3(8, 9, 10).toString()).toEqual('(8.00; 9.00; 10.00)');
    expect(new Vector3(8, 9, 10).toString(0)).toEqual('(8; 9; 10)');
})

test('check toArray property', () => {
    expect(new Vector3(8, 9, 10).toArray()).toEqual([8, 9, 10]);
})

