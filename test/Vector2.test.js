import Vector2 from "../src/Vector2";

test('check getting coordinates from vector', () => {
    let vec0 = Vector2.zero;
    let vec1 = new Vector2(2, 3);

    expect(vec0.x).toBe(0);
    expect(vec0.y).toBe(0);
    expect(vec1.x).toBe(2);
    expect(vec1.y).toBe(3);
})

test('default vectors', () => {
    expect(Vector2.zero.x).toEqual(0);
    expect(Vector2.zero.y).toEqual(0);

    expect(Vector2.right.x).toEqual(1);
    expect(Vector2.right.y).toEqual(0);

    expect(Vector2.up.x).toEqual(0);
    expect(Vector2.up.y).toEqual(1);
})

test('check magnitude and magnitude squared calculation of the vector', () => {
    let vec1 = new Vector2(1, 1);
    let vec2 = new Vector2(0, 0);
    let vec3 = new Vector2(-2, -3);
    let vec4 = new Vector2(Math.sqrt(9), Math.sqrt(7));

    expect(vec1.magnitude).toEqual(Math.sqrt(2));
    expect(vec1.magnitudeSquared).toEqual(2);
    expect(vec2.magnitude).toEqual(0);
    expect(vec2.magnitudeSquared).toEqual(0);
    expect(vec3.magnitude).toEqual(Math.sqrt(13));
    expect(vec3.magnitudeSquared).toEqual(13);
    expect(vec4.magnitude).toEqual(4);
    expect(vec4.magnitudeSquared).toEqual(16);
})

test('check normalize property', () => {
    let vec1 = new Vector2(1, 1);
    let vec3 = new Vector2(-2, -3);
    let vec4 = new Vector2(7, 8);

    expect(vec1.normalize().x).toEqual(1 / Math.sqrt(2));
    expect(vec1.normalize().y).toEqual(1 / Math.sqrt(2));

    expect(vec3.normalize().x).toEqual(vec3.x / vec3.magnitude);
    expect(vec3.normalize().y).toEqual(vec3.y / vec3.magnitude);

    expect(vec4.normalize().x).toEqual(vec4.x / vec4.magnitude);
    expect(vec4.normalize().y).toEqual(vec4.y / vec4.magnitude);
})

test('check neg property', () => {
    expect(new Vector2(9, -2).neg().x).toEqual(-9);
    expect(new Vector2(9, -2).neg().y).toEqual(2);

    expect(new Vector2(0, 0).neg().x).toEqual(-0);
})

test('check add function', () => {
    let vec0 = Vector2.zero;
    let vec1 = new Vector2(-2, 1);
    let vec2 = new Vector2(-4, 5);
    let vecSum = Vector2.add(new Vector2(1, -1), new Vector2(-5, 4));
    expect(Vector2.add(vec0, vec0)).toEqual(vec0);

    expect(Vector2.add(vec1, vec2).x).toEqual(-6);
    expect(Vector2.add(vec1, vec2).y).toEqual(6);

    expect(Vector2.add(vecSum, vec1).x).toEqual(-6);
    expect(Vector2.add(vecSum, vec1).y).toEqual(4);
})

test('check sub function', () => {
    let vec0 = Vector2.zero;
    let vec1 = new Vector2(-5, -9);
    let vec2 = new Vector2(-4, 5);
    let vecSum = Vector2.add(vec1, vec2);
    expect(Vector2.sub(vec1, vec1)).toEqual(Vector2.zero);
    expect(Vector2.sub(vec0, vec2)).toEqual(new Vector2(4, -5));
    expect(Vector2.sub(vecSum, vec1)).toEqual(vec2);
})

test('check mul function', () => {
    let vec0 = Vector2.zero;
    let vec1 = new Vector2(4, 3);
    let vec2 = new Vector2(-2, -3);
    let vec3 = new Vector2(-4, 5);
    expect(Vector2.mul(vec0, Math.random())).toEqual(Vector2.zero);
    expect(Vector2.mul(vec1, 0)).toEqual(Vector2.zero);
    expect(Vector2.mul(vec1, -1)).toEqual(new Vector2(-4, -3));
    expect(Vector2.mul(vec2, -1)).toEqual(new Vector2(2, 3));
    expect(Vector2.mul(vec3, 3)).toEqual(new Vector2(-12, 15));
    expect(Vector2.mul(vec2, -10)).toEqual(new Vector2(20, 30));
})

test('check div function', () => {
    let vec0 = Vector2.zero;
    let vec1 = new Vector2(4, 3);
    let vec2 = new Vector2(-2, -3);
    let vec3 = new Vector2(-4, 5);

    expect(Vector2.div(vec0, 0)).toEqual(new Vector2(null, null));
    expect(Vector2.div(vec1, 0)).toEqual(new Vector2(null, null));
    expect(Vector2.div(vec0, 2)).toEqual(vec0);
    expect(Vector2.div(vec2, -2)).toEqual(new Vector2(1, 1.5));
    expect(Vector2.div(vec3, -20)).toEqual(new Vector2(0.2, -0.25));
    expect(Vector2.div(vec3, 20)).toEqual(new Vector2(-0.2, 0.25));
})

test('check dot product function', () => {
    let vec0 = new Vector2(0, 2);

    expect(Vector2.dot(vec0, new Vector2(1, 0))).toEqual(0);
    expect(Vector2.dot(new Vector2(-6, -7), new Vector2(2, -3))).toEqual(9);
    expect(Vector2.dot(vec0, new Vector2(9, 1))).toEqual(2);
})

test('check cross product function', () => {
    let vec0 = new Vector2(0, 2);

    expect(Vector2.cross(vec0, new Vector2(0, 1))).toEqual(0);
    expect(Vector2.cross(new Vector2(-6, -7), new Vector2(2, -3))).toEqual(32);
    expect(Vector2.cross(vec0, new Vector2(9, 1))).toEqual(-18);
})

test('check angleBetween function', () => {
    let vec0 = new Vector2(0, 2);
    let vec1 = new Vector2(1, 0);
    let vec2 = new Vector2(2, -3);
    let vec3 = new Vector2(4, -6);
    let vec4 = new Vector2(1, -8);

    expect(Vector2.angleBetween(vec0, vec1)).toEqual(Math.PI / 2);
    expect(Vector2.angleBetween(vec2, vec3)).toEqual(0);
    expect(Vector2.angleBetween(vec2, vec4)).toBeCloseTo(
        Math.acos(26 / (Math.sqrt(13) * Math.sqrt(65))),
        1);
    expect(Vector2.angleBetween(vec3, vec4)).toBeCloseTo(
        Math.acos(52 / (Math.sqrt(52) * Math.sqrt(65))),
        1);
})

test('check distance function', () => {
    let vec0 = new Vector2(1, 2);
    let vec1 = new Vector2(-4, 0);
    let vec2 = new Vector2(9, -8);

    expect(typeof Vector2.distance(vec0, vec1)).toEqual('number');
    expect(Vector2.distance(vec0, vec0)).toEqual(0);
    expect(Vector2.distance(vec0, vec1)).toEqual(Math.sqrt(29));
    expect(Vector2.distance(vec0, vec2)).toEqual(Math.sqrt(164));
    expect(Vector2.distance(vec1, vec2)).toEqual(Math.sqrt(233));
})

test('check linearly interpolation between two vectors (lerp function)', () => {
    let vec0 = new Vector2(-2, 4);
    let vec1 = new Vector2(1, -5);

    expect(Vector2.lerp(vec0, vec1, 0.6)).toEqual(new Vector2(-0.2, -1.4));
    expect(Vector2.lerp(vec0, vec1, -0.6)).toEqual(new Vector2(-2, 4));
    expect(Vector2.lerp(vec0, vec1, 2)).toEqual(new Vector2(1, -5));
    expect(Vector2.lerp(new Vector2(0, 1), new Vector2(1, 0), 0.5)).toEqual(new Vector2(0.5, 0.5));
})

test('check fromArray method', () => {
    expect(Vector2.fromArray([0, 2]) instanceof Vector2).toBeTruthy();
    expect(Vector2.fromArray([0, 2])).toEqual(new Vector2(0, 2));
})

test('check toString property', () => {
    expect(typeof new Vector2(8, 9).toString()).toEqual('string');
    expect(new Vector2(8, 9).toString()).toEqual('(8.00; 9.00)');
    expect(new Vector2(8, 9).toString(0)).toEqual('(8; 9)');
})

test('check toArray property', () => {
    expect(new Vector2(8, 9).toArray()).toEqual([8, 9]);
})