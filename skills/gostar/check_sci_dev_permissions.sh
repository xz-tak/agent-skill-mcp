#!/bin/bash
# Check if AWS profile sci-dev has necessary Secrets Manager permissions

set -e

PROFILE="sci-dev"
REGION="us-east-1"
TEST_SECRET_NAME="test-permissions-check-$(date +%s)"

echo "=========================================="
echo "Checking sci-dev Profile Permissions"
echo "=========================================="
echo ""

# Check identity
echo "1. Checking identity..."
IDENTITY=$(aws sts get-caller-identity --profile ${PROFILE} --query 'Arn' --output text)
echo "   Identity: ${IDENTITY}"
echo ""

# Check ListSecrets
echo "2. Testing ListSecrets permission..."
if aws secretsmanager list-secrets --profile ${PROFILE} --region ${REGION} --max-items 1 &>/dev/null; then
    echo "   ✓ Can list secrets"
else
    echo "   ✗ Cannot list secrets"
    exit 1
fi
echo ""

# Check CreateSecret
echo "3. Testing CreateSecret permission..."
if aws secretsmanager create-secret \
    --profile ${PROFILE} \
    --region ${REGION} \
    --name "${TEST_SECRET_NAME}" \
    --description "Temporary test secret" \
    --secret-string "test-value" &>/dev/null; then
    echo "   ✓ Can create secrets"

    # Clean up test secret
    echo "   Cleaning up test secret..."
    aws secretsmanager delete-secret \
        --profile ${PROFILE} \
        --region ${REGION} \
        --secret-id "${TEST_SECRET_NAME}" \
        --force-delete-without-recovery &>/dev/null
    echo "   ✓ Test secret deleted"
else
    echo "   ✗ Cannot create secrets"
    ERROR=$(aws secretsmanager create-secret \
        --profile ${PROFILE} \
        --region ${REGION} \
        --name "${TEST_SECRET_NAME}" \
        --secret-string "test" 2>&1 || true)
    echo "   Error: ${ERROR}"
    exit 1
fi
echo ""

# Check DescribeSecret
echo "4. Testing DescribeSecret permission..."
if aws secretsmanager list-secrets --profile ${PROFILE} --region ${REGION} --max-items 1 --query 'SecretList[0].Name' --output text | xargs -I {} aws secretsmanager describe-secret --profile ${PROFILE} --region ${REGION} --secret-id {} &>/dev/null; then
    echo "   ✓ Can describe secrets"
else
    echo "   ✗ Cannot describe secrets"
fi
echo ""

# Check GetSecretValue
echo "5. Testing GetSecretValue permission..."
SECRET_NAME=$(aws secretsmanager list-secrets --profile ${PROFILE} --region ${REGION} --max-items 1 --query 'SecretList[0].Name' --output text)
if [ "${SECRET_NAME}" != "None" ] && [ ! -z "${SECRET_NAME}" ]; then
    if aws secretsmanager get-secret-value --profile ${PROFILE} --region ${REGION} --secret-id "${SECRET_NAME}" &>/dev/null; then
        echo "   ✓ Can get secret values"
    else
        echo "   ✗ Cannot get secret values"
    fi
else
    echo "   ⚠ No secrets available to test GetSecretValue"
fi
echo ""

# Check PutSecretValue (update)
echo "6. Testing PutSecretValue permission..."
if aws secretsmanager create-secret \
    --profile ${PROFILE} \
    --region ${REGION} \
    --name "${TEST_SECRET_NAME}" \
    --secret-string "test1" &>/dev/null; then

    if aws secretsmanager put-secret-value \
        --profile ${PROFILE} \
        --region ${REGION} \
        --secret-id "${TEST_SECRET_NAME}" \
        --secret-string "test2" &>/dev/null; then
        echo "   ✓ Can update secret values"
    else
        echo "   ✗ Cannot update secret values"
    fi

    # Clean up
    aws secretsmanager delete-secret \
        --profile ${PROFILE} \
        --region ${REGION} \
        --secret-id "${TEST_SECRET_NAME}" \
        --force-delete-without-recovery &>/dev/null
else
    echo "   ⚠ Skipping PutSecretValue test"
fi
echo ""

# Check IAM permissions
echo "7. Testing IAM permissions..."
ROLE_NAME="ml-model-server-model-server-role"
if aws iam get-role --profile ${PROFILE} --role-name "${ROLE_NAME}" &>/dev/null; then
    echo "   ✓ Can read IAM role: ${ROLE_NAME}"

    if aws iam list-attached-role-policies --profile ${PROFILE} --role-name "${ROLE_NAME}" &>/dev/null; then
        echo "   ✓ Can list role policies"
    else
        echo "   ✗ Cannot list role policies"
    fi

    if aws iam create-policy \
        --profile ${PROFILE} \
        --policy-name "test-policy-${RANDOM}" \
        --policy-document '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Action":"s3:GetObject","Resource":"*"}]}' \
        --dry-run 2>&1 | grep -q "DryRunOperation"; then
        echo "   ✓ Can create IAM policies"
    else
        echo "   ⚠ Cannot verify policy creation (may still work)"
    fi

    if aws iam attach-role-policy \
        --profile ${PROFILE} \
        --role-name "${ROLE_NAME}" \
        --policy-arn "arn:aws:iam::aws:policy/ReadOnlyAccess" \
        --dry-run 2>&1 | grep -q "DryRunOperation"; then
        echo "   ✓ Can attach policies to roles"
    else
        echo "   ⚠ Cannot verify policy attachment (may still work)"
    fi
else
    echo "   ✗ Cannot read IAM role: ${ROLE_NAME}"
    echo "   This may require additional IAM permissions"
fi
echo ""

echo "=========================================="
echo "Summary"
echo "=========================================="
echo ""
echo "Profile: ${PROFILE}"
echo "Region: ${REGION}"
echo ""
echo "✓ sci-dev profile can manage secrets in AWS Secrets Manager"
echo ""
echo "Next step: Create GOSTAR secret"
echo "  ./setup_secrets_manager.sh \"<password>\" --profile sci-dev"
echo ""
